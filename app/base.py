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
import tempfile
import sys
import traceback

app = Flask(__name__, template_folder='templates', static_folder='static')

# Ensure Jinja2 is set up correctly
app.jinja_env.auto_reload = True
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.secret_key = 'your_secret_key_here'  # Required for session

# Directory to store saved results
SAVE_DIR = '/var/www/html/students_25/Team7/app/saved_results'
try:
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR)
    test_file = os.path.join(SAVE_DIR, 'test_write.txt')
    with open(test_file, 'w') as f:
        f.write('test')
    print(f"Successfully verified write permissions to {SAVE_DIR}")
except Exception as e:
    print(f"Warning: Could not write to {SAVE_DIR}: {str(e)}")
    SAVE_DIR = tempfile.mkdtemp(prefix='genomic_results_')
    print(f"Using temporary directory instead: {SAVE_DIR}")
    
SAVE_METADATA = os.path.join(SAVE_DIR, 'saved_files.json')
saved_files = {}

def connect_database(hostname='bioed-new.bu.edu', port=4253, database='Team7', 
                    username='', password=''):
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

def execute_query(cursor, condition_name, cell_type, gene_params, 
                  output_fields, cre_fields, tf_fields, include_de=False, 
                  de_params=None, cre_params=None, tf_params=None,
                  page=1, per_page=50):
    """Execute queries based on parameters with pagination."""
    # Build the base query
    query_parts = ["SELECT DISTINCT"]
    select_fields = []
    
    # Add requested gene fields
    if len(output_fields) > 0:
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
                select_fields.append("bp.name as pathway")
    
    # Add DE fields if requested
    if include_de and de_params:
        de_fields = de_params.get('de_fields', [])
        for field in de_fields:
            select_fields.append(f"de.{field}")
    
    # Add related CRE info if requested
    if len(cre_fields) > 0:
        for field in cre_fields:
            if field == 'cre_chr':
                select_fields.append("cre.chromosome as cre_chr")
            elif field == 'cre_start':
                select_fields.append("cre.start_position as cre_start")
            elif field == 'cre_end':
                select_fields.append("cre.end_position as cre_end")
            elif field == 'cre_log2fc':
                select_fields.append("cre.cre_log2foldchange as cre_log2fc")
            elif field == 'cre_padj':
                select_fields.append("cre.padj as cre_padj")
            elif field == 'cre_distance':
                select_fields.append("cgi.distance_to_TSS as cre_distance")
    
    if len(tf_fields) > 0:
        for field in tf_fields:
            if field == 'tf_checkbox':
                select_fields.append("tf.name as tf")            

                
    # Base cases
    if not select_fields and len(output_fields) > 0: 
        select_fields = ["g.gene_symbol as hgnc_symbol", "g.Entrez_ID as entrez_id", "g.chromosome", "g.start_position", "g.end_position"]
    elif not select_fields and len(output_fields) == 0 and len(cre_fields) > 0:
        select_fields = ["cre.chromosome as cre_chr", "cre.start_position as cre_start", "cre.end_position as cre_end", "cre.cre_log2foldchange as cre_log2fc"]
    elif not select_fields and len(output_fields) == 0 and len(cre_fields) == 0 and len(tf_fields) > 0:
        select_fields = ["tf.name as tf"]
    elif not select_fields and len(output_fields) == 0 and len(cre_fields) == 0 and len(tf_fields) == 0:
        select_fields = ["c.name as condition_name", "ct.cell as cell_type"]
    
    query_parts.append(", ".join(select_fields))   
    
    # Basic gene query 
    query_parts.append("""
    FROM Genes g
    JOIN Differential_Expression de ON g.gid = de.gid
    JOIN Conditions c ON de.cdid = c.cdid AND c.name = %s
    JOIN Cell_Type ct ON de.cell_id = ct.cell_id AND ct.cell = %s
    JOIN Gene_Pathway_Associations gpa ON g.gid = gpa.gid
    JOIN Biological_Pathways bp ON gpa.pid = bp.pid
    """)
    
    # Add related CRE info 
    if len(cre_fields) > 0 or len(cre_params) > 0:
        query_parts.append("""
        JOIN CRE_Gene_Interactions cgi ON g.gid = cgi.gid
        JOIN Cis_Regulatory_Elements cre ON cgi.cid = cre.cid AND cre.cdid = c.cdid AND cre.cell_id = ct.cell_id
        """)
        
    # safeguard if user bypasses cre to get tfs
    if (len(tf_fields)+len(tf_params)) > 0 and (len(cre_fields)+len(cre_params)) == 0:
        query_parts.append("""
        JOIN CRE_Gene_Interactions cgi ON g.gid = cgi.gid
        JOIN Cis_Regulatory_Elements cre ON cgi.cid = cre.cid AND cre.cdid = c.cdid AND cre.cell_id = ct.cell_id
        """)
        if len(cre_fields) == 0:
            cre_fields = ["cre_chr", "cre_start", "cre_end"]
        if len(cre_params) == 0:
            cre_params = {}
        
    # include tfs
    if len(tf_fields) > 0 or len(tf_params) > 0:
        query_parts.append("""
        JOIN Merged_CRES mc ON cre.mcid = mc.mcid
        JOIN TF_CRE_Interactions tci ON mc.mcid = tci.mcid AND tci.cdid = c.cdid AND tci.cell_id = ct.cell_id
        JOIN Transcription_Factors tf ON tci.tfid = tf.tfid
        """)
    
    
    params = [condition_name, cell_type]
    # Add WHERE clause
    query_parts.append("WHERE 1=1")
    
    # Add gene-specific filters
    if len(gene_params) > 0:
        if gene_params.get('gene-id-type') and gene_params.get('gene-identifier'):
            id_type = gene_params.get('gene-id-type')
            identifier = gene_params.get('gene-identifier')
            if id_type == 'hgnc':
                query_parts.append("AND lower(g.gene_symbol) = lower(%s)")
                params.append(identifier)
            elif id_type == 'entrez':
                query_parts.append("AND g.Entrez_ID = %s")
                params.append(identifier)
            elif id_type == 'ensembl':
                query_parts.append("AND g.Ensembl_ID = %s")
                params.append(identifier)

    # Add chromosome, position, and pathway filters if provided
    if len(gene_params) > 0:
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
            pathway_input = gene_params.get('gene-pathway').upper()
            query_parts.append("AND lower(bp.name) LIKE lower(%s)")
            params.append(f"%{pathway_input}%")

    # Add DE filters if requested
    if include_de and de_params:
        if de_params.get('padj_filter'):
            query_parts.append("AND de.padj < %s")
            params.append(float(de_params.get('padj_filter')))
        
        if de_params.get('logfc_filter'):
            query_parts.append("AND abs(de.log2foldchange) > %s")
            params.append(float(de_params.get('logfc_filter')))
            
    if len(cre_params) > 0:
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
            
    if len(tf_params) > 0:
        # Add TF-specific filters
        if tf_params.get('tf-name'):
            query_parts.append("AND lower(tf.name) = lower(%s)")
            params.append(tf_params.get('tf-name'))
    
    # Combine all parts into a base query (without pagination)
    base_query = " ".join(query_parts)
    
    # First get total count for pagination metadata
    count_query = f"SELECT COUNT(*) FROM ({base_query}) as count_query"
    try:
        cursor.execute(count_query, params)
        total_count = cursor.fetchone()[0]
    except mariadb.Error as e:
        return None, None, f"Database count query error: {str(e)}"
    
    # Add pagination to the final query
    paginated_query = base_query + f" LIMIT {per_page} OFFSET {(page-1)*per_page}"
    
    # Execute the paginated query
    try:
        cursor.execute(paginated_query, params)
        results = cursor.fetchall()
        # Convert results to list of dictionaries with column names
        column_names = [desc[0] for desc in cursor.description]
        result_dicts = [dict(zip(column_names, row)) for row in results]
        
        # Create pagination metadata
        pagination_info = {
            'total_records': total_count,
            'page': page,
            'per_page': per_page,
            'total_pages': (total_count + per_page - 1) // per_page  # Ceiling division
        }
        
        return result_dicts, pagination_info, None
    except mariadb.Error as e:
        return None, None, f"Database query error: {str(e)}\nQuery: {paginated_query}\nParams: {params}"

def generate_table_html(results, headers, pagination_info, title=None, **params):
    """Generate HTML table from results with server-side pagination."""
    if not results:
        return "<p>No results found matching your criteria.</p>"
    
    # Extract pagination information
    page = pagination_info['page']
    per_page = pagination_info['per_page']  # Should be 10
    total_records = pagination_info['total_records']
    total_pages = pagination_info['total_pages']

    # Generate hidden fields for pagination
    hidden_fields = []
    for key, value in params.items():
        if key != 'page' and key != 'per_page':  
            if isinstance(value, dict):  
                for sub_key, sub_value in value.items():
                    if sub_value is not None:
                        hidden_fields.append(f'<input type="hidden" name="{key}[{sub_key}]" value="{sub_value}">')
            elif isinstance(value, list):  
                for item in value:
                    hidden_fields.append(f'<input type="hidden" name="{key}[]" value="{item}">')
            elif value is not None:
                hidden_fields.append(f'<input type="hidden" name="{key}" value="{value}">')
    
    # Add hidden fields for specific parameter types that require special handling
    if 'output_fields' in params and params['output_fields']:
        for field in params['output_fields']:
            hidden_fields.append(f'<input type="hidden" name="output-fields" value="{field}">')
    
    if 'cre_fields' in params and params['cre_fields']:
        for field in params['cre_fields']:
            hidden_fields.append(f'<input type="hidden" name="cre-output-fields" value="{field}">')
    
    if 'tf_fields' in params and params['tf_fields']:
        for field in params['tf_fields']:
            hidden_fields.append(f'<input type="hidden" name="tf-checkbox" value="{field}">')
    
    if 'include_de' in params and params['include_de']:
        hidden_fields.append('<input type="hidden" name="include_de" value="on">')
    
    if 'de_params' in params and params['de_params']:
        de_params = params['de_params']
        if 'de_fields' in de_params and de_params['de_fields']:
            for field in de_params['de_fields']:
                hidden_fields.append(f'<input type="hidden" name="de_fields" value="{field}">')
        
        if 'padj_filter' in de_params and de_params['padj_filter']:
            hidden_fields.append(f'<input type="hidden" name="padj_filter" value="{de_params["padj_filter"]}">')
        
        if 'logfc_filter' in de_params and de_params['logfc_filter']:
            hidden_fields.append(f'<input type="hidden" name="logfc_filter" value="{de_params["logfc_filter"]}">')
    
    if 'gene_params' in params and params['gene_params']:
        gene_params = params['gene_params']
        for key, value in gene_params.items():
            if value:
                hidden_fields.append(f'<input type="hidden" name="{key}" value="{value}">')
    
    if 'cre_params' in params and params['cre_params']:
        cre_params = params['cre_params']
        for key, value in cre_params.items():
            if value:
                hidden_fields.append(f'<input type="hidden" name="{key}" value="{value}">')
    
    if 'tf_params' in params and params['tf_params']:
        tf_params = params['tf_params']
        for key, value in tf_params.items():
            if value:
                hidden_fields.append(f'<input type="hidden" name="{key}" value="{value}">')
    
    # Always include condition and cell_type in hidden fields
    if 'condition' in params and params['condition']:
        hidden_fields.append(f'<input type="hidden" name="condition" value="{params["condition"]}">')
    
    if 'cell_type' in params and params['cell_type']:
        hidden_fields.append(f'<input type="hidden" name="cell_type" value="{params["cell_type"]}">')
    
    if 'active_tab' in params and params['active_tab']:
        hidden_fields.append(f'<input type="hidden" name="active_tab" value="{params["active_tab"]}">')
    
    # Create the results table
    table_html = f"""
    <div class="results-section">
        <h2>{title or "Search Results"}</h2>
        <div class="results-meta">
            <span class="badge results-count">{len(results)} results</span>
            <span class="badge page-info">Page {page} of {total_pages}</span>
        </div>
        
        <div class="table-container">
            <table class="results-table">
                <thead>
                    <tr>
    """
    
    for header in headers:
        table_html += f"<th>{header.replace('_', ' ').title()}</th>"
    
    table_html += """
            </tr></thead>
            <tbody>
    """
    
    for row in results:
        table_html += "<tr>"
        for header in headers:
            value = row.get(header, '')
            table_html += f"<td>{value if value is not None else ''}</td>"
        table_html += "</tr>"
    
    table_html += """
            </tbody>
        </table>
    </div>
    """

    # Add server-side pagination controls
    search_url = url_for('search')
    
    # Only show pagination if there are multiple pages
    if total_pages > 1:
        pagination = f"""
        <div class="pagination">
            <a href="javascript:void(0)" 
            class="pagination-link {'disabled' if page <= 1 else ''}" 
            data-page="{max(1, page-1)}">Previous</a>
            <span>Page {page} of {total_pages}</span>
            <a href="javascript:void(0)" 
            class="pagination-link {'disabled' if page >= total_pages else ''}" 
            data-page="{min(total_pages, page+1)}">Next</a>

            <span id="loading-spinner" style="display:none; margin-left:10px;">
                <i class="fas fa-spinner fa-spin"></i>
            </span>
        </div>
        """
        return table_html + pagination
    else:
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
    
    try:
        # Write to CSV
        with open(filepath, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(results)
        return True
    except Exception as e:
        print(f"Error saving CSV: {str(e)}")
        return False

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

def load_saved_files():
    """Load saved files metadata from JSON file."""
    try:
        if os.path.exists(SAVE_METADATA):
            with open(SAVE_METADATA, 'r') as f:
                return json.load(f)
    except Exception as e:
        print(f"Error loading saved files metadata: {str(e)}")
    return {}

def save_saved_files(saved_files_dict):
    """Save files metadata to JSON file."""
    try:
        with open(SAVE_METADATA, 'w') as f:
            json.dump(saved_files_dict, f, indent=2)
        return True
    except Exception as e:
        print(f"Error saving files metadata: {str(e)}")
        return False

# Modified route: Change the index route to render base.html instead
@app.route('/')
def index():
    # Render base.html as the home page
    return render_template('base.html')

# Add a dedicated route for the search page
@app.route('/search_page', methods=['GET'])
def search_page():
    return render_template('updated_search.html', 
                         table_html=None,
                         condition=None,
                         cell_type=None,
                         active_tab='gene',
                         error=None,
                         result_id=None)

# Add routes for all other HTML pages
@app.route('/guide')
def guide():
    return render_template('guide.html')

@app.route('/faq')
def faq():
    return render_template('faq.html')

@app.route('/resources')
def resources():
    return render_template('Resources_DownloadData.html')

@app.route('/aboutus')
def aboutus():
    return render_template('about_us.html')

@app.route('/visualizations')
def visualizations():
    return render_template('visualizations.html')

@app.route('/github')
def github():
    return render_template('github.html')

@app.route('/contactus')
def contactus():
    return render_template('contact_us.html')

@app.route('/citation')
def citation():
    return render_template('citation.html')

@app.route('/search', methods=['GET'])
def search():
    # Debug information
    print(f"Request Method: {request.method}")
    print(f"Args Data: {request.args}")
    
    # Check if it's a search request (has parameters)
    if not request.args:
        return redirect(url_for('search_page'))
    
    # Get data from the request object
    data = request.args
    
    # Get common parameters
    condition = data.get('condition')
    cell_type = data.get('cell_type')
    active_tab = data.get('active_tab', 'gene')  # Default to gene tab
    
    # Get pagination parameters
    try:
        page = int(data.get('page', 1))
    except (ValueError, TypeError):
        page = 1
        
    try:
        per_page = int(data.get('per_page', 10))  # Default to 10 items per page
    except (ValueError, TypeError):
        per_page = 10
    
    # Print debug info
    print(f"Pagination: page={page}, per_page={per_page}")
    print(f"Params: condition={condition}, cell_type={cell_type}, active_tab={active_tab}")
    
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
    
    # new!! add for ajax
    is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
    
    try:
        results = None
        error = None
        table_html = ""
        
        # Generate a unique ID for this result set
        result_id = str(uuid.uuid4())
        search_type = active_tab
        save_option = data.get('save_option', 'view')  # Options: view, save, both
        
        if active_tab == 'gene' or active_tab == 'cre' or active_tab == 'tf':
            # Get gene parameters
            gene_params = {
                'gene-id-type': data.get('gene-id-type'),
                'gene-identifier': data.get('gene-identifier'),
                'gene-chr': data.get('gene-chr'),
                'gene-start': data.get('gene-start'),
                'gene-end': data.get('gene-end'),
                'gene-pathway': data.get('gene-pathway')
            }
            
            # Get output fields from GET request
            output_fields = data.getlist('output-fields')
            
            # Get differential expression parameters
            include_de = data.get('include_de') == 'on'
            de_params = None
            if include_de:
                de_params = {
                    'de_fields': data.getlist('de_fields'),
                    'padj_filter': data.get('padj_filter'),
                    'logfc_filter': data.get('logfc_filter')
                }
            
            # Get CRE parameters
            cre_params = {
                'cre-chr': data.get('cre-chr'),
                'cre-start': data.get('cre-start'),
                'cre-end': data.get('cre-end'),
                'cre-log2fc': data.get('cre-log2fc')
            }
            
            # Get CRE fields from GET request
            cre_fields = data.getlist('cre-output-fields')
            
            # Get TF parameters
            tf_params = {
                'tf-name': data.get('tf-name')  
            }
            
            # Get TF fields from GET request
            tf_fields = data.getlist('tf-checkbox')
            
            # Call execute_query with pagination parameters
            results, pagination_info, error = execute_query(
                cursor, condition, cell_type, gene_params,
                output_fields, cre_fields, tf_fields, include_de=include_de, 
                de_params=de_params, cre_params=cre_params, tf_params=tf_params,
                page=page, per_page=per_page
            )
            
            # Build the title for results
            title = f"Search Results ({condition}, {cell_type})"
            description = f"Search for {condition} in {cell_type} cells"
            if gene_params.get('gene-identifier'):
                title += f" - {gene_params.get('gene-id-type').upper()}: {gene_params.get('gene-identifier')}"
                description += f", {gene_params.get('gene-id-type').upper()}: {gene_params.get('gene-identifier')}"
            
            # Generate column headers and table HTML for the result table
            if results:
                headers = list(results[0].keys())
                table_html = generate_table_html(
                        results=results, 
                        headers=headers, 
                        pagination_info=pagination_info, 
                        title=title,
                        condition=condition,
                        cell_type=cell_type,
                        active_tab=active_tab,
                        gene_params=gene_params,
                        output_fields=output_fields,
                        include_de=include_de,
                        de_params=de_params,
                        cre_params=cre_params,
                        cre_fields=cre_fields,
                        tf_params=tf_params,
                        tf_fields=tf_fields
                    )

                if is_ajax:
                    return render_template('results_fragment.html',
                                            table_html=table_html,
                                            results=results,
                                            headers=headers,
                                            pagination_info=pagination_info,
                                            title=title,
                                            condition=condition,
                                            cell_type=cell_type,
                                            active_tab=active_tab,
                                            result_id=result_id if results else None)
                                            # gene_params=gene_params,
                                            # output_fields=output_fields,
                                            # include_de=include_de,
                                            # de_params=de_params,
                                            # cre_params=cre_params,
                                            # cre_fields=cre_fields,
                                            # tf_params=tf_params,
                                            # tf_fields=tf_fields,
                                            # error=error)
                else:
                    # original thml rendering for full page loads
                    table_html = table_html
                    return render_template('updated_search.html',
                                           table_html=table_html,
                                           condition=condition,
                                           cell_type=cell_type,
                                           active_tab=active_tab,
                                           error=error,
                                           result_id=result_id if results else None,
                                           pagination_info=pagination_info)
            else:
                if is_ajax:
                    return jsonify({
                        'status': 'error',
                        'message': error or "No results found matching your criteria."
                    })
                else:
                    return render_template('updated_search.html',
                                           error=error or "No results found matching your criteria.",
                                           table_html=None,
                                           condition=None,
                                           cell_type=None,
                                           active_tab=active_tab,
                                           result_id=None)
    except Exception as e:
        if is_ajax:
            return jsonify({
                'status': 'error',
                'message': f"Application error: {str(e)}"
            }), 500
        else:
            return render_template('updated_search.html',
                                  error=f"Application error: {str(e)}",
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
    # Load saved files metadata
    saved_files = load_saved_files()
    
    # Convert saved_files dict to list for template
    files_list = list(saved_files.values())
    return render_template('downloads.html', saved_files=files_list)

@app.route('/download/<file_id>')
def download_file(file_id):
    """Download a specific saved file."""
    # Load saved files metadata
    saved_files = load_saved_files()
    
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
    # Load saved files metadata
    saved_files = load_saved_files()
    
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
    
    # Save updated metadata
    save_saved_files(saved_files)
    
    return "File deleted", 200

@app.route('/save_current_result/<result_id>', methods=['POST'])
def save_current_result(result_id):
    """Save currently displayed results to downloads."""
    try:
        # Check if results exist in the session
        if 'current_results' not in session:
            print("No results found in session")
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
            
            title = f"{search_type.capitalize()} Results ({condition}, {cell_type})"
            description = f"Search for {condition} in {cell_type} cells"
            
            # Load current saved files
            saved_files = load_saved_files()
            
            # Add the new file
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
            
            # Save updated metadata
            save_saved_files(saved_files)
            
            return redirect(url_for('downloads'))
        else:
            print("Failed to save results to CSV")
            return "Failed to save results", 500
    except Exception as e:
        print(f"Error in save_current_result: {str(e)}")
        import traceback
        traceback.print_exc()
        return f"Error: {str(e)}", 500

@app.route('/students_25/yhkwok/HW3_folder/yhkwok_visualization/volcano_plot', methods=['POST'])
@app.route('/volcano_plot', methods=['POST'])
def volcano_plot():
    if request.method == 'POST':
        condition_name = request.form.get('condition_name')
        cell_type = request.form.get('cell_type')
        
        if not condition_name or not cell_type:
            return jsonify([])
        
        conn = None
        cursor = None
        try:
            conn = connect_database()
            cursor = conn.cursor(dictionary=True)  
            
            query = """
            SELECT g.gene_symbol, de.log2foldchange, de.p_value, de.padj
            FROM Differential_Expression de
            JOIN Genes g ON de.gid = g.gid
            JOIN Conditions c ON de.cdid = c.cdid
            JOIN Cell_Type ct ON de.cell_id = ct.cell_id
            WHERE c.name = ? AND ct.cell = ?
            ORDER BY g.gene_symbol
            """
            
            cursor.execute(query, (condition_name, cell_type))
            results = cursor.fetchall()
            
            return jsonify(results)
        
        except Exception as e:
            return jsonify({"error": f"Database error occurred: {str(e)}"}), 500
            
        finally:
            if cursor:
                cursor.close()
            if conn:
                conn.close()
            
    return jsonify([])

@app.route('/students_25/yhkwok/HW3_folder/yhkwok_visualization/fgsea_plot', methods=['POST'])
@app.route('/fgsea_plot', methods=['POST'])
def fgsea_plot():
    if request.method == 'POST':
        condition_name = request.form.get('condition_name')
        cell_type = request.form.get('cell_type')
        pathway_count = request.form.get('pathway_count', 10)
        
        try:
            pathway_count = int(pathway_count)
            if pathway_count < 1:
                pathway_count = 10
            elif pathway_count > 50:
                pathway_count = 50
        except ValueError:
            pathway_count = 10
            
        if not condition_name or not cell_type:
            return jsonify([])
        
        conn = None
        cursor = None
        try:
            conn = connect_database()
            cursor = conn.cursor(dictionary=True)
            
            up_query = """
            SELECT 
                bp.name AS pathway_name, 
                COUNT(DISTINCT g.gid) AS gene_count,
                SUM(CASE WHEN de.log2foldchange > 0 THEN 1 ELSE 0 END) AS up_regulated,
                SUM(CASE WHEN de.log2foldchange < 0 THEN 1 ELSE 0 END) AS down_regulated,
                AVG(de.log2foldchange) AS avg_fold_change,
                -LOG10(GREATEST(MIN(COALESCE(de.padj, 1)), 0.000001)) AS neg_log_padj,
                'up' AS regulation_direction
            FROM Biological_Pathways bp
            JOIN Gene_Pathway_Associations gpa ON bp.pid = gpa.pid
            JOIN Genes g ON gpa.gid = g.gid
            JOIN Differential_Expression de ON g.gid = de.gid
            JOIN Conditions c ON de.cdid = c.cdid
            JOIN Cell_Type ct ON de.cell_id = ct.cell_id
            WHERE 
                c.name = ? 
                AND ct.cell = ?
                AND de.padj < 0.05
            GROUP BY bp.name
            HAVING 
                COUNT(DISTINCT g.gid) >= 3
                AND SUM(CASE WHEN de.log2foldchange > 0 THEN 1 ELSE 0 END) > SUM(CASE WHEN de.log2foldchange < 0 THEN 1 ELSE 0 END)
            ORDER BY neg_log_padj DESC
            LIMIT ?
            """
            
            down_query = """
            SELECT 
                bp.name AS pathway_name, 
                COUNT(DISTINCT g.gid) AS gene_count,
                SUM(CASE WHEN de.log2foldchange > 0 THEN 1 ELSE 0 END) AS up_regulated,
                SUM(CASE WHEN de.log2foldchange < 0 THEN 1 ELSE 0 END) AS down_regulated,
                AVG(de.log2foldchange) AS avg_fold_change,
                -LOG10(GREATEST(MIN(COALESCE(de.padj, 1)), 0.000001)) AS neg_log_padj,
                'down' AS regulation_direction
            FROM Biological_Pathways bp
            JOIN Gene_Pathway_Associations gpa ON bp.pid = gpa.pid
            JOIN Genes g ON gpa.gid = g.gid
            JOIN Differential_Expression de ON g.gid = de.gid
            JOIN Conditions c ON de.cdid = c.cdid
            JOIN Cell_Type ct ON de.cell_id = ct.cell_id
            WHERE 
                c.name = ? 
                AND ct.cell = ?
                AND de.padj < 0.05
            GROUP BY bp.name
            HAVING 
                COUNT(DISTINCT g.gid) >= 3
                AND SUM(CASE WHEN de.log2foldchange > 0 THEN 1 ELSE 0 END) <= SUM(CASE WHEN de.log2foldchange < 0 THEN 1 ELSE 0 END)
            ORDER BY neg_log_padj DESC
            LIMIT ?
            """
            
            cursor.execute(up_query, (condition_name, cell_type, pathway_count))
            up_results = cursor.fetchall()
            
            cursor.execute(down_query, (condition_name, cell_type, pathway_count))
            down_results = cursor.fetchall()
            
            results = up_results + down_results
            
            return jsonify(results)
        
        except Exception as e:
            return jsonify({"error": f"Database error occurred: {str(e)}"}), 500
            
        finally:
            if cursor:
                cursor.close()
            if conn:
                conn.close()
            
    return jsonify([])

@app.route('/students_25/yhkwok/HW3_folder/yhkwok_visualization/cre_gene_scatter', methods=['POST'])
@app.route('/cre_gene_scatter', methods=['POST'])
def cre_gene_scatter():
    if request.method == 'POST':
        condition_name = request.form.get('condition_name')
        cell_type = request.form.get('cell_type')
        
        if not condition_name or not cell_type:
            return jsonify([])
        
        conn = None
        cursor = None
        try:
            conn = connect_database()
            cursor = conn.cursor(dictionary=True)
            
            query = """
            SELECT 
                g.gene_symbol, 
                de.log2foldchange as gene_log2fc, 
                cre.cre_log2foldchange as cre_log2fc,
                de.padj as gene_padj,
                0.05 as cre_padj,
                cgi.distance_to_TSS,
                cre.chromosome as cre_chr,
                cre.start_position as cre_start,
                cre.end_position as cre_end
            FROM Genes g
            JOIN Differential_Expression de ON g.gid = de.gid
            JOIN Conditions c ON de.cdid = c.cdid AND c.name = ?
            JOIN Cell_Type ct ON de.cell_id = ct.cell_id AND ct.cell = ?
            JOIN CRE_Gene_Interactions cgi ON g.gid = cgi.gid
            JOIN Cis_Regulatory_Elements cre ON cgi.cid = cre.cid
                AND cre.cdid = c.cdid 
                AND cre.cell_id = ct.cell_id
            ORDER BY g.gene_symbol
            """
            
            cursor.execute(query, (condition_name, cell_type))
            results = cursor.fetchall()
            
            return jsonify(results)
        
        except Exception as e:
            try:
                fallback_query = """
                SELECT 
                    g.gene_symbol, 
                    de.log2foldchange as gene_log2fc, 
                    cre.cre_log2foldchange as cre_log2fc,
                    de.padj as gene_padj,
                    0.05 as cre_padj,
                    cgi.distance_to_TSS,
                    cre.chromosome as cre_chr,
                    cre.start_position as cre_start,
                    cre.end_position as cre_end
                FROM Genes g
                JOIN Differential_Expression de ON g.gid = de.gid
                JOIN Conditions c ON de.cdid = c.cdid 
                JOIN Cell_Type ct ON de.cell_id = ct.cell_id
                JOIN CRE_Gene_Interactions cgi ON g.gid = cgi.gid
                JOIN Cis_Regulatory_Elements cre ON cgi.cid = cre.cid
                WHERE c.name = ? AND ct.cell = ?
                ORDER BY g.gene_symbol
                LIMIT 100
                """
                cursor.execute(fallback_query, (condition_name, cell_type))
                results = cursor.fetchall()
                return jsonify(results)
            except Exception:
                return jsonify({"error": f"Database error occurred: {str(e)}"}), 500
            
        finally:
            if cursor:
                cursor.close()
            if conn:
                conn.close()
            
    return jsonify([])

@app.route('/students_25/yhkwok/HW3_folder/yhkwok_visualization/get_conditions', methods=['GET'])
@app.route('/get_conditions', methods=['GET'])
def get_conditions():
    conn = None
    cursor = None
    try:
        conn = connect_database()
        cursor = conn.cursor(dictionary=True)
        
        query = "SELECT name FROM Conditions ORDER BY name"
        cursor.execute(query)
        results = cursor.fetchall()
        
        return jsonify([item['name'] for item in results])
    except Exception as e:
        return jsonify({"error": f"Database error occurred: {str(e)}"}), 500
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()
            
@app.route('/students_25/yhkwok/HW3_folder/yhkwok_visualization/get_cell_types', methods=['GET'])
@app.route('/get_cell_types', methods=['GET'])
def get_cell_types():
    conn = None
    cursor = None
    try:
        conn = connect_database()
        cursor = conn.cursor(dictionary=True)
        
        query = "SELECT cell FROM Cell_Type ORDER BY cell"
        cursor.execute(query)
        results = cursor.fetchall()
        
        return jsonify([item['cell'] for item in results])
    except Exception as e:
        return jsonify({"error": f"Database error occurred: {str(e)}"}), 500
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()
        
@app.route('/test_db_connection', methods=['GET'])
def test_db_connection():
    try:
        conn = connect_database()
        cursor = conn.cursor()
        cursor.execute("SELECT 1")
        result = cursor.fetchone()
        cursor.close()
        conn.close()
        return jsonify({"status": "success", "message": "Database connection successful", "result": result})
    except Exception as e:
        return jsonify({"status": "error", "message": f"Database connection failed: {str(e)}"}), 500


if __name__ == '__main__':
    app.run(debug=True)
