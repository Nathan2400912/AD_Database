#!/usr/bin/env python3

from flask import Flask, request, render_template, jsonify, redirect, url_for, send_file, make_response, session
import mariadb
from string import Template
import json
import os
import csv
import io
import datetime
import uuid
import shutil
from werkzeug.utils import secure_filename

app = Flask(__name__, template_folder='templates', static_folder='static')

# Ensure Jinja2 is set up correctly
app.jinja_env.auto_reload = True
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.secret_key = 'your_secret_key_here'  # Required for session

# Directory to store saved results
SAVE_DIR = 'saved_results'
if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

# In-memory storage for file metadata
# In production, this would be in a database
saved_files = {}

def connect_database(hostname='localhost', port=3306, database='geneinteract', 
                    username='geneuser', password='password'):
    """Connect to the MariaDB database."""
    try:
        connection = mariadb.connect(
            host=hostname,
            user=username,
            password=password,
            db=database,
            port=int(port)
        )
        cursor = connection.cursor()
        return connection, cursor
    except mariadb.Error as e:
        return None, str(e)

def execute_gene_query(cursor, condition_name, cell_type, gene_params, output_fields, include_de=False, de_params=None, include_related=False):
    """Execute gene-related query based on parameters."""
    # Build the base query
    query_parts = ["SELECT"]
    select_fields = []
    
    # Add requested gene fields
    for field in output_fields:
        if field == 'hgnc':
            select_fields.append("g.gene_symbol as hgnc_symbol")
        elif field == 'entrez':
            select_fields.append("g.Entrez_ID as entrez_id")
        elif field == 'ensembl':
            select_fields.append("g.Ensembl_ID as ensembl_id")
        elif field == 'chr':
            select_fields.append("g.chromosome")
        elif field == 'start':
            select_fields.append("g.start_position")
        elif field == 'end':
            select_fields.append("g.end_position")
        elif field == 'strand':
            select_fields.append("g.strand")
        elif field == 'pathway':
            select_fields.append("GROUP_CONCAT(DISTINCT bp.name SEPARATOR ', ') as pathway")
    
    # If no fields selected, select default fields
    if not select_fields:
        select_fields = ["g.gene_symbol as hgnc_symbol", "g.Entrez_ID as entrez_id", "g.chromosome", "g.start_position", "g.end_position"]
    
    # Add DE fields if requested
    if include_de and de_params:
        de_fields = de_params.get('de_fields', [])
        for field in de_fields:
            select_fields.append(f"de.{field}")
    
    # Add related CRE info if requested
    if include_related:
        select_fields.append("GROUP_CONCAT(DISTINCT CONCAT(cre.chromosome, ':', cre.start_position, '-', cre.end_position) SEPARATOR '; ') as associated_cres")
        select_fields.append("GROUP_CONCAT(DISTINCT tf.name SEPARATOR ', ') as associated_tfs")
    
    query_parts.append(", ".join(select_fields))
    
    # Get condition and cell_type IDs
    query_parts.append("""
    FROM Genes g
    JOIN Conditions c ON c.name = %s
    JOIN Cell_Type ct ON ct.cell = %s
    """)
    
    params = [condition_name, cell_type]
    
    # Add DE join if requested
    if include_de or include_related:
        query_parts.append("""
        LEFT JOIN Differential_Expression de ON g.gid = de.gid AND de.cdid = c.cdid AND de.cell_id = ct.cell_id
        """)
    
    # Add related CRE and TF info if requested
    if include_related:
        query_parts.append("""
        LEFT JOIN CRE_Gene_Interactions cgi ON g.gid = cgi.gid
        LEFT JOIN Cis_Regulatory_Elements cre ON cgi.cid = cre.cid AND cre.cdid = c.cdid AND cre.cell_id = ct.cell_id
        LEFT JOIN TF_CRE_Interactions tci ON cre.mcid = tci.mcid AND tci.cdid = c.cdid AND tci.cell_id = ct.cell_id
        LEFT JOIN Transcription_Factors tf ON tci.tfid = tf.tfid
        """)
    
    # Add pathway join
    query_parts.append("""
    LEFT JOIN Gene_Pathway_Associations gpa ON g.gid = gpa.gid
    LEFT JOIN Biological_Pathways bp ON gpa.pid = bp.pid
    """)
    
    # Add WHERE clause
    query_parts.append("WHERE 1=1")
    
    # Add gene-specific filters
    if gene_params.get('gene-id-type') and gene_params.get('gene-identifier'):
        id_type = gene_params.get('gene-id-type')
        identifier = gene_params.get('gene-identifier')
        if id_type == 'hgnc':
            query_parts.append("AND g.gene_symbol = %s")
            params.append(identifier)
        elif id_type == 'entrez':
            query_parts.append("AND g.Entrez_ID = %s")
            params.append(identifier)
        elif id_type == 'ensembl':
            query_parts.append("AND g.Ensembl_ID = %s")
            params.append(identifier)
    
    # Add chromosome, position, and pathway filters if provided
    if gene_params.get('gene-chr'):
        query_parts.append("AND g.chromosome = %s")
        params.append(gene_params.get('gene-chr'))
    
    if gene_params.get('gene-start'):
        query_parts.append("AND g.start_position >= %s")
        params.append(int(gene_params.get('gene-start')))
    
    if gene_params.get('gene-end'):
        query_parts.append("AND g.end_position <= %s")
        params.append(int(gene_params.get('gene-end')))
    
    if gene_params.get('gene-pathway'):
        query_parts.append("AND bp.name LIKE %s")
        params.append(f"%{gene_params.get('gene-pathway')}%")
    
    # Add DE filters if requested
    if include_de and de_params:
        if de_params.get('padj_filter'):
            query_parts.append("AND de.padj < %s")
            params.append(float(de_params.get('padj_filter')))
        
        if de_params.get('logfc_filter'):
            query_parts.append("AND abs(de.log2foldchange) > %s")
            params.append(float(de_params.get('logfc_filter')))
    
    # Add GROUP BY clause for proper aggregation
    query_parts.append("GROUP BY g.gid")
    
    # Combine all parts into a single query
    query = " ".join(query_parts)
    
    # Execute the query
    try:
        cursor.execute(query, params)
        results = cursor.fetchall()
        # Convert results to list of dictionaries with column names
        column_names = [desc[0] for desc in cursor.description]
        result_dicts = [dict(zip(column_names, row)) for row in results]
        return result_dicts, None
    except mariadb.Error as e:
        return None, f"Database query error: {str(e)}\nQuery: {query}\nParams: {params}"

def execute_cre_query(cursor, condition_name, cell_type, cre_params, output_fields, include_related=False):
    """Execute CRE-related query based on parameters."""
    # Build the base query
    query_parts = ["SELECT"]
    select_fields = []
    
    # Add requested CRE fields
    for field in output_fields:
        if field == 'chr':
            select_fields.append("cre.chromosome")
        elif field == 'start':
            select_fields.append("cre.start_position")
        elif field == 'end':
            select_fields.append("cre.end_position")
        elif field == 'log2fc':
            select_fields.append("cre.cre_log2foldchange")
        elif field == 'padj':
            select_fields.append("cre.padj")
        elif field == 'activity':
            select_fields.append("cre.activity_score")
    
    # If no fields selected, select default fields
    if not select_fields:
        select_fields = ["cre.chromosome", "cre.start_position", "cre.end_position", "cre.cre_log2foldchange"]
    
    # Add related gene and TF info if requested
    if include_related:
        select_fields.append("GROUP_CONCAT(DISTINCT g.gene_symbol SEPARATOR ', ') as associated_genes")
        select_fields.append("GROUP_CONCAT(DISTINCT tf.name SEPARATOR ', ') as associated_tfs")
        select_fields.append("GROUP_CONCAT(DISTINCT CONCAT(g.gene_symbol, ' (', cgi.distance_to_TSS, 'bp)') SEPARATOR '; ') as gene_distances")
    
    query_parts.append(", ".join(select_fields))
    
    # Get condition and cell_type IDs
    query_parts.append("""
    FROM Cis_Regulatory_Elements cre
    JOIN Conditions c ON c.name = %s
    JOIN Cell_Type ct ON ct.cell = %s
    """)
    
    params = [condition_name, cell_type]
    
    # Add related gene and TF info if requested
    if include_related:
        query_parts.append("""
        LEFT JOIN CRE_Gene_Interactions cgi ON cre.cid = cgi.cid
        LEFT JOIN Genes g ON cgi.gid = g.gid
        LEFT JOIN TF_CRE_Interactions tci ON cre.mcid = tci.mcid AND tci.cdid = c.cdid AND tci.cell_id = ct.cell_id
        LEFT JOIN Transcription_Factors tf ON tci.tfid = tf.tfid
        """)
    
    # Add WHERE clause for condition and cell_type
    query_parts.append("WHERE cre.cdid = c.cdid AND cre.cell_id = ct.cell_id")
    
    # Add CRE-specific filters
    if cre_params.get('cre-chr'):
        query_parts.append("AND cre.chromosome = %s")
        params.append(cre_params.get('cre-chr'))
    
    if cre_params.get('cre-start'):
        query_parts.append("AND cre.start_position >= %s")
        params.append(int(cre_params.get('cre-start')))
    
    if cre_params.get('cre-end'):
        query_parts.append("AND cre.end_position <= %s")
        params.append(int(cre_params.get('cre-end')))
    
    if cre_params.get('cre-log2fc'):
        query_parts.append("AND abs(cre.cre_log2foldchange) > %s")
        params.append(float(cre_params.get('cre-log2fc')))
    
    # Add GROUP BY clause for proper aggregation if using related data
    if include_related:
        query_parts.append("GROUP BY cre.cid")
    
    # Combine all parts into a single query
    query = " ".join(query_parts)
    
    # Execute the query
    try:
        cursor.execute(query, params)
        results = cursor.fetchall()
        # Convert results to list of dictionaries with column names
        column_names = [desc[0] for desc in cursor.description]
        result_dicts = [dict(zip(column_names, row)) for row in results]
        return result_dicts, None
    except mariadb.Error as e:
        return None, f"Database query error: {str(e)}\nQuery: {query}\nParams: {params}"

def execute_tf_query(cursor, condition_name, cell_type, tf_params, include_related=False):
    """Execute TF-related query based on parameters."""
    # Build the base query
    query_parts = ["SELECT tf.name as tf_name"]
    
    # Add count of binding sites and related info if requested
    if include_related:
        query_parts.append(", COUNT(DISTINCT tci.mcid) as binding_sites_count")
        query_parts.append(", GROUP_CONCAT(DISTINCT CONCAT(cre.chromosome, ':', cre.start_position, '-', cre.end_position) SEPARATOR '; ') as bound_cres")
        query_parts.append(", GROUP_CONCAT(DISTINCT g.gene_symbol SEPARATOR ', ') as regulated_genes")
    
    # Join tables
    query_parts.append("""
    FROM Transcription_Factors tf
    JOIN Conditions c ON c.name = %s
    JOIN Cell_Type ct ON ct.cell = %s
    JOIN TF_CRE_Interactions tci ON tf.tfid = tci.tfid AND tci.cdid = c.cdid AND tci.cell_id = ct.cell_id
    """)
    
    params = [condition_name, cell_type]
    
    # Add related CRE and gene info if requested
    if include_related:
        query_parts.append("""
        JOIN Cis_Regulatory_Elements cre ON tci.mcid = cre.mcid AND cre.cdid = c.cdid AND cre.cell_id = ct.cell_id
        LEFT JOIN CRE_Gene_Interactions cgi ON cre.cid = cgi.cid
        LEFT JOIN Genes g ON cgi.gid = g.gid
        """)
    
    # Add TF name filter if provided
    if tf_params.get('tf-name'):
        query_parts.append("WHERE tf.name = %s")
        params.append(tf_params.get('tf-name'))
    else:
        query_parts.append("WHERE 1=1")
    
    # Add GROUP BY clause
    query_parts.append("GROUP BY tf.tfid")
    
    # Add ORDER BY binding sites count if related data is included
    if include_related:
        query_parts.append("ORDER BY binding_sites_count DESC")
    
    # Combine all parts into a single query
    query = " ".join(query_parts)
    
    # Execute the query
    try:
        cursor.execute(query, params)
        results = cursor.fetchall()
        # Convert results to list of dictionaries with column names
        column_names = [desc[0] for desc in cursor.description]
        result_dicts = [dict(zip(column_names, row)) for row in results]
        return result_dicts, None
    except mariadb.Error as e:
        return None, f"Database query error: {str(e)}\nQuery: {query}\nParams: {params}"

def generate_table_html(results, headers, title=None):
    """Generate HTML table from results."""
    if not results:
        return "<p>No results found matching your criteria.</p>"
    
    # Create HTML table
    table_html = """
    <div class="results-section">
        <h2>{title}</h2>
        <p>Found {count} records matching your search criteria.</p>
        <div class="table-container">
            <table>
                <thead>
                    <tr>
    """.format(title=title or "Search Results", count=len(results))
    
    # Add table headers
    for header in headers:
        # Format header for display (convert snake_case or camelCase to Title Case)
        display_header = header.replace('_', ' ').title()
        table_html += f"<th>{display_header}</th>"
    
    table_html += """
                    </tr>
                </thead>
                <tbody>
    """
    
    # Add table rows
    for row in results:
        table_html += "<tr>"
        for header in headers:
            cell_value = row.get(header, '')
            if cell_value is None:
                cell_value = ''
            table_html += f"<td>{cell_value}</td>"
        table_html += "</tr>"
    
    table_html += """
                </tbody>
            </table>
        </div>
    </div>
    """
    
    return table_html

def save_results_to_csv(results, filename):
    """Save query results to a CSV file."""
    if not results:
        return False
    
    # Create directory if it doesn't exist
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    
    # Get the headers from the first result
    headers = list(results[0].keys())
    
    filepath = os.path.join(SAVE_DIR, filename)
    
    # Write to CSV
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(results)
    
    return True

def get_file_size(filepath):
    """Get human-readable file size."""
    size_bytes = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"

def get_preview(filepath, max_rows=5):
    """Get a preview of the CSV file."""
    preview_rows = []
    with open(filepath, 'r') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # Get headers
        for i, row in enumerate(reader):
            if i >= max_rows:
                break
            preview_rows.append(row)
    
    # Truncate long text
    for i, row in enumerate(preview_rows):
        for j, cell in enumerate(row):
            if len(cell) > 50:
                preview_rows[i][j] = cell[:50] + "..."
    
    if not preview_rows:
        return "No data available for preview."
    
    return ", ".join([f"{headers[0]}: {row[0]}" for row in preview_rows])

# Modified route: Change the index route to render base.html instead
@app.route('/', methods=['GET'])
def index():
    # Render base.html as the home page
    return render_template('base.html')

# Add a dedicated route for the search page
@app.route('/search_page', methods=['GET'])
def search_page():
    # This will render the search page without executing a search
    return render_template('updated_search.html', 
                         table_html=None,
                         condition=None,
                         cell_type=None,
                         active_tab='gene',
                         error=None,
                         result_id=None)

# Add routes for all HTML pages
@app.route('/base')
def base():
    return render_template('base.html')

@app.route('/guide')
def guide():
    return render_template('guide.html')

@app.route('/faq')
def faq():
    return render_template('faq.html')

@app.route('/resources')
def resources():
    return render_template('Resources_DownloadData.html')

@app.route('/search', methods=['GET', 'POST'])
def search():
    # Check if it's a search request (has parameters)
    if request.method == 'GET' and not request.args:
        return redirect(url_for('search_page'))
    
    data = request.form if request.method == 'POST' else request.args
    
    # Get common parameters
    condition = data.get('condition')
    cell_type = data.get('cell_type')
    active_tab = data.get('active_tab', 'gene')  # Default to gene tab
    
    # Validate required parameters
    if not condition or not cell_type:
        error = "Error: Both condition and cell type are required."
        return render_template('updated_search.html', 
                              error=error,
                              table_html=None,
                              condition=None,
                              cell_type=None,
                              active_tab=active_tab,
                              result_id=None)
    
    #Connect to database
    connection, cursor = connect_database()
    if not connection:
        error_message = f"Error: Could not connect to the database. {cursor}"
        return render_template('updated_search.html', 
                              error=error_message,
                              table_html=None,
                              condition=None,
                              cell_type=None,
                              active_tab=active_tab,
                              result_id=None)
    
    try:
        results = None
        error = None
        table_html = ""
        
        # Generate a unique ID for this result set
        result_id = str(uuid.uuid4())
        search_type = active_tab
        save_option = data.get('save_option', 'view')  # Options: view, save, both
        
        # Execute query based on active tab
        if active_tab == 'gene':
            # Parse gene parameters
            gene_params = {
                'gene-id-type': data.get('gene-id-type'),
                'gene-identifier': data.get('gene-identifier'),
                'gene-chr': data.get('gene-chr'),
                'gene-start': data.get('gene-start'),
                'gene-end': data.get('gene-end'),
                'gene-pathway': data.get('gene-pathway')
            }
            
            # Parse output fields
            output_fields = request.form.getlist('output-fields') if request.method == 'POST' else request.args.getlist('output-fields')
            
            # Check if related data is requested
            include_related = data.get('include-related-gene') == 'on'
            
            # Check if DE data is requested
            include_de = data.get('include_de') == 'on'
            de_params = None
            if include_de:
                de_params = {
                    'de_fields': request.form.getlist('de_fields') if request.method == 'POST' else request.args.getlist('de_fields'),
                    'padj_filter': data.get('padj_filter'),
                    'logfc_filter': data.get('logfc_filter')
                }
            
            # Execute gene query
            results, error = execute_gene_query(cursor, condition, cell_type, gene_params, output_fields, include_de, de_params, include_related)
            
            # Build the title for results
            title = f"Gene Search Results ({condition}, {cell_type})"
            description = f"Gene search for {condition} in {cell_type} cells"
            if gene_params.get('gene-identifier'):
                title += f" - {gene_params.get('gene-id-type').upper()}: {gene_params.get('gene-identifier')}"
                description += f", {gene_params.get('gene-id-type').upper()}: {gene_params.get('gene-identifier')}"
            
            # Generate column headers for the result table
            if results:
                headers = list(results[0].keys())
                table_html = generate_table_html(results, headers, "Gene Search Results")
            else:
                table_html = "<p>No gene results found matching your criteria.</p>"
        
        elif active_tab == 'cre':
            # Parse CRE parameters
            cre_params = {
                'cre-chr': data.get('cre-chr'),
                'cre-start': data.get('cre-start'),
                'cre-end': data.get('cre-end'),
                'cre-log2fc': data.get('cre-log2fc')
            }
            
            # Parse output fields
            output_fields = request.form.getlist('cre-output-fields') if request.method == 'POST' else request.args.getlist('cre-output-fields')
            
            # Check if related data is requested
            include_related = data.get('include-related-cre') == 'on'
            
            # Execute CRE query
            results, error = execute_cre_query(cursor, condition, cell_type, cre_params, output_fields, include_related)
            
            # Build the title for results
            title = f"CRE Search Results ({condition}, {cell_type})"
            description = f"Cis-regulatory element search for {condition} in {cell_type} cells"
            if cre_params.get('cre-chr'):
                chrom = cre_params.get('cre-chr')
                start = cre_params.get('cre-start', 'start')
                end = cre_params.get('cre-end', 'end')
                title += f" - {chrom}:{start}-{end}"
                description += f", region {chrom}:{start}-{end}"
            
            # Generate column headers for the result table
            if results:
                headers = list(results[0].keys())
                table_html = generate_table_html(results, headers, "CRE Search Results")
            else:
                table_html = "<p>No CRE results found matching your criteria.</p>"
        
        elif active_tab == 'tf':
            # Parse TF parameters
            tf_params = {
                'tf-name': data.get('tf-name')
            }
            
            # Check if related data is requested
            include_related = data.get('include-related-tf') == 'on'
            
            # Execute TF query
            results, error = execute_tf_query(cursor, condition, cell_type, tf_params, include_related)
            
            # Build the title for results
            title = f"TF Search Results ({condition}, {cell_type})"
            description = f"Transcription factor search for {condition} in {cell_type} cells"
            if tf_params.get('tf-name'):
                title += f" - {tf_params.get('tf-name')}"
                description += f", factor: {tf_params.get('tf-name')}"
            
            # Generate column headers for the result table
            if results:
                headers = list(results[0].keys())
                table_html = generate_table_html(results, headers, "Transcription Factor Search Results")
            else:
                table_html = "<p>No transcription factor results found matching your criteria.</p>"
        
        else:
            error = "Error: Invalid search type."
            table_html = ""
        
        # Save the results to a file if requested
        if results and (save_option == 'save' or save_option == 'both'):
            current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{search_type}_{condition}_{cell_type}_{current_time}_{result_id}.csv"
            
            if save_results_to_csv(results, filename):
                # Store file metadata
                filepath = os.path.join(SAVE_DIR, filename)
                file_size = get_file_size(filepath)
                preview = get_preview(filepath)
                
                saved_files[result_id] = {
                    'id': result_id,
                    'filename': filename,
                    'filepath': filepath,
                    'title': title,
                    'description': description,
                    'type': search_type.capitalize(),
                    'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    'size': file_size,
                    'preview': preview,
                    'condition': condition,
                    'cell_type': cell_type
                }
                
                # If save only, redirect to downloads page
                if save_option == 'save':
                    return redirect(url_for('downloads'))
            else:
                error = "Error: Failed to save results to file."
        
        # Store results in session for save_current_result
        if results:
            session['current_results'] = results
        
        # Return results for display
        return render_template('updated_search.html', 
                              table_html=table_html,
                              condition=condition,
                              cell_type=cell_type,
                              active_tab=active_tab,
                              error=error,
                              result_id=result_id if results else None)
        
    except Exception as e:
        error_message = f"Application error: {str(e)}"
        return render_template('updated_search.html', 
                              error=error_message,
                              table_html=None,
                              condition=None,
                              cell_type=None,
                              active_tab=active_tab,
                              result_id=None)
        
    finally:
        if connection:
            cursor.close()
            connection.close()

@app.route('/downloads')
def downloads():
    """Display the downloads page with saved files."""
    # Convert saved_files dict to list for template
    files_list = list(saved_files.values())
    return render_template('downloads.html', saved_files=files_list)

@app.route('/download/<file_id>')
def download_file(file_id):
    """Download a specific saved file."""
    if file_id not in saved_files:
        return "File not found", 404
    
    file_data = saved_files[file_id]
    return send_file(file_data['filepath'], 
                    mimetype='text/csv',
                    as_attachment=True,
                    download_name=file_data['filename'])

@app.route('/delete-file/<file_id>', methods=['DELETE'])
def delete_file(file_id):
    """Delete a saved file."""
    if file_id not in saved_files:
        return "File not found", 404
    
    # Remove file from filesystem
    filepath = saved_files[file_id]['filepath']
    try:
        os.remove(filepath)
    except OSError:
        pass  # File may not exist
    
    # Remove file from metadata
    del saved_files[file_id]
    
    return "File deleted", 200

@app.route('/combine-files', methods=['POST'])
def combine_files():
    """Combine multiple files into one CSV."""
    data = request.json
    file_ids = data.get('fileIds', [])
    
    if not file_ids or len(file_ids) < 2:
        return "At least two files must be selected", 400
    
    # Check if all files exist
    for file_id in file_ids:
        if file_id not in saved_files:
            return f"File {file_id} not found", 404
    
    # Create a new combined file
    combined_id = str(uuid.uuid4())
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"combined_{current_time}_{combined_id}.csv"
    filepath = os.path.join(SAVE_DIR, filename)
    
    # Combine files
    combined_data = []
    headers_set = set()
    
    # First pass: collect all headers
    for file_id in file_ids:
        file_path = saved_files[file_id]['filepath']
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for header in reader.fieldnames:
                headers_set.add(header)
    
    headers = sorted(list(headers_set))
    
    # Second pass: combine data
    for file_id in file_ids:
        file_path = saved_files[file_id]['filepath']
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Create a new row with all headers
                new_row = {header: row.get(header, '') for header in headers}
                combined_data.append(new_row)
    
    # Write combined data
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(combined_data)
    
    # Store combined file metadata
    description = "Combined file from "
    description += ", ".join([saved_files[file_id]['title'] for file_id in file_ids[:3]])
    if len(file_ids) > 3:
        description += f" and {len(file_ids) - 3} more"
    
    file_size = get_file_size(filepath)
    preview = get_preview(filepath)
    
    saved_files[combined_id] = {
        'id': combined_id,
        'filename': filename,
        'filepath': filepath,
        'title': f"Combined Results ({len(combined_data)} records)",
        'description': description,
        'type': 'Combined',
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'size': file_size,
        'preview': preview,
        'condition': 'Multiple',
        'cell_type': 'Multiple'
    }
    
    # Return the combined file for download
    return send_file(filepath, 
                    mimetype='text/csv',
                    as_attachment=True,
                    download_name=filename)

@app.route('/save-current/<result_id>', methods=['POST'])
def save_current_result(result_id):
    """Save currently displayed results to downloads."""
    # Check if results exist in the session
    if not hasattr(request, 'session') or 'current_results' not in session:
        return "No results to save", 400
    
    results = session['current_results']
    
    # Generate filename and save
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    search_type = request.form.get('search_type', 'query')
    condition = request.form.get('condition', 'unknown')
    cell_type = request.form.get('cell_type', 'unknown')
    
    filename = f"{search_type}_{condition}_{cell_type}_{current_time}_{result_id}.csv"
    
    if save_results_to_csv(results, filename):
        # Store file metadata
        filepath = os.path.join(SAVE_DIR, filename)
        file_size = get_file_size(filepath)
        preview = get_preview(filepath)
        
        title = request.form.get('title', f"{search_type.capitalize()} Results")
        description = request.form.get('description', f"Search for {condition} in {cell_type} cells")
        
        saved_files[result_id] = {
            'id': result_id,
            'filename': filename,
            'filepath': filepath,
            'title': title,
            'description': description,
            'type': search_type.capitalize(),
            'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'size': file_size,
            'preview': preview,
            'condition': condition,
            'cell_type': cell_type
        }
        
        return redirect(url_for('downloads'))
    else:
        return "Failed to save results", 500

if __name__ == '__main__':
    app.run(debug=True)