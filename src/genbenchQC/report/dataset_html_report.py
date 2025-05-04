from datetime import datetime

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HTML Report Output</title>
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-3.0.1.min.js" charset="utf-8"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: flex;
            max-width: 100%; /* Prevents content from exceeding the viewport width */
        }
        .sidebar {
            width: 250px;
            background: #f4f4f4;
            padding: 20px;
            box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
            height: 100vh; /* Full height */
            position: fixed; /* Fixed position */
            overflow-y: auto;
            z-index: 1000; /* Ensures it stays above other elements */
        }
        .sidebar a {
            display: block;
            margin: 10px 0;
            text-decoration: none;
            color: #333;
        }
        .content {
            margin-left: 300px;
            padding: 20px;
            width: calc(100% - 350px); /* Adjust width to avoid overlap */
            overflow-x: hidden;
        }
        section {
            margin-bottom: 50px;
            width: 100%;
        }
        .grid-container {
            display: grid;
            grid-template-columns: 1fr 1fr;
            grid-gap: 20px;
        }
        #sequence-duplication-levels table {
            table-layout: fixed; /* Ensures columns respect defined widths */
            width: 100%; /* Makes the table take up the full width of its container */
            border-collapse: collapse; /* Optional: Makes the table look cleaner */
        }

        #sequence-duplication-levels table th,
        #sequence-duplication-levels table td {
            padding: 8px; /* Adds padding for better readability */
            border: 1px solid #ddd; /* Adds a border for clarity */
        }

        #sequence-duplication-levels table td {
            white-space: nowrap; /* Prevents text from wrapping */
            overflow: hidden; /* Hides overflowing text */
            text-overflow: ellipsis; /* Adds "..." to indicate clipped text */
        }

        #sequence-duplication-levels table th {
            text-align: left; /* Aligns header text to the left */
        }

        .count_column {
            width: 5%; /* Fixed width for the first column */
        }

        .sequence_column {
            width: 95%; /* Fixed width for the second column */
        }

        h1 {
            text-align: center;
            margin-bottom: 50px;
        }

        h2 {
            color: #333;
            border-bottom: 2px solid #ddd;
            padding-bottom: 5px;
        }

        h3 {
            color: #555;
            margin-bottom: 15px;
        }

        .chart-container {
            margin-bottom: 30px;
        }

        .data-item {
            font-size: 1.2em;
            margin-bottom: 10px;
        }

        .data-item span {
            font-weight: bold;
            font-size: 1em;
        }
    </style>
</head>
<body>

<!-- DEFINE MAIN PAGE LAYOUT -->

    <div class="container">

<!-- SIDEBAR - NAVIGATION -->
<!-- This sidebar will be fixed on the left side of the page -->

        <div class="sidebar">
            <h2>Navigation</h2>
            <a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a>
            <div style="margin-left: 15px;"> <!-- Using 15px margin for indentation -->
                <a href="#filename">Filename</a>
                <a href="#num-sequences">Number of sequences</a>
                <a href="#num-bases">Number of bases</a>
                <a href="#unique-bases">Unique bases</a>
                <a href="#gc-content">%GC content</a>
                <a href="#dedup-sequences">Number of sequences left after deduplication</a>
            </div>
            <a href="#general-descriptive-statistics">General Descriptive Statistics</a>
            <div style="margin-left: 15px;"> <!-- Using 15px margin for indentation -->
                <a href="#sequence-lengths">Sequence lengths</a>
                <a href="#duplicit-sequences">Duplicit sequences</a>
            </div>
            <a href="#per-sequence-descriptive-stats">Per Sequence Descriptive Stats</a>
            <div style="margin-left: 15px;"> <!-- Using 15px margin for indentation -->
                <a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a>
                <a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a>
                <a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a>
                <a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a>
                <a href="#per-sequence-gc-content">Per Sequence GC Content</a>
            </div>
        </div> <!-- END OF SIDEBAR -->

    
<!-- MAIN CONTENT AREA -->
<!-- This is where the main content of the page will be displayed -->
<!-- Using section tags to define different sections of the report -->
<!-- Using divs to create placeholders for the data and charts -->

        <div class="content">
            <h1>HTML Report Output</h1>

            <section id="basic-descriptive-statistics">
                <h2>Basic Descriptive Statistics</h2>
                <div class="grid-container">
                    <div class="left-column">
                        <div class="data-item" id="d1-filename">
                            <span>Filename: </span> 
                            {{d1-filename}} <!-- Filename will be displayed here -->
                        </div>
                        <div class="data-item" id="d1-num-sequences">
                            <span>Number of sequences: </span> 
                            {{d1-number_of_sequences}} <!-- Number of sequences will be displayed here -->
                        </div>
                        <div class="data-item" id="d1-num-bases">
                            <span>Number of bases: </span> 
                            {{d1-number_of_bases}} <!-- Number of bases will be displayed here -->
                        </div>
                        <div class="data-item" id="d1-unique-bases">
                            <span>Unique bases: </span> 
                            {{d1-unique_bases}} <!-- Unique bases will be displayed here -->
                        </div>
                        <div class="data-item" id="d1-gc-content">
                            <span>%GC content: </span> 
                            {{d1-gc_content}} <!-- %GC content will be displayed here -->
                        </div>
                        <div class="data-item" id="d1-dedup-sequences">
                            <span>Number of sequences left after deduplication:</span> 
                            {{d1-dedup_sequences}} <!-- Number of sequences left after deduplication will be displayed here -->
                        </div>
                    </div>
                    <div class="right-column">
                        <div class="data-item" id="d2-filename">
                            <span>Filename:</span> 
                            {{d2-filename}}
                        </div>
                        <div class="data-item" id="d2-num-sequences">
                            <span>Number of sequences:</span> 
                            {{d2-number_of_sequences}}
                        </div>
                        <div class="data-item" id="d2-num-bases">
                            <span>Number of bases:</span> 
                            {{d2-number_of_bases}}
                        </div>
                        <div class="data-item" id="d2-unique-bases">
                            <span>Unique bases:</span> 
                            {{d2-unique_bases}}
                        </div>
                        <div class="data-item" id="d2-gc-content">
                            <span>%GC content: </span>
                            {{d2-gc_content}}
                        </div>
                        <div class="data-item" id="d2-dedup-sequences">
                            <span>Number of sequences left after deduplication:</span> 
                            {{d2-dedup_sequences}}
                        </div>
                    </div>
                </div> <!-- END OF GRID CONTAINER -->
            </section> <!-- END OF BASIC DESCRIPTIVE STATISTICS SECTION -->

            <section id="general-descriptive-statistics">
                <h2>General Descriptive Statistics</h2>

                <h3>Sequence lengths</h3>
                <div class="chart-container" id="sequence-lengths">

                </div>

                <h3>Duplicit sequences</h3>

                <div class="chart-container" id="duplicit-sequences">
                    
                </div>
            </section> <!-- END OF GENERAL DESCRIPTIVE STATISTICS SECTION -->

            <section id="per-sequence-descriptive-stats">
                <h2>Per Sequence Descriptive Stats</h2>

                <h3>Per Sequence Nucleotide Content</h3>
                <div class="chart-container" id="per-sequence-nucleotide-content">
                </div>

                <h3>Per Sequence Dinucleotide Content</h3>
                <div class="chart-container" id="per-sequence-dinucleotide-content">
                </div>

                <h3>Per Position Nucleotide Content</h3>
                <div class="chart-container" id="per-position-nucleotide-content">
                </div>

                <h3>Per Position Reversed Nucleotide Content</h3>
                <div class="chart-container" id="per-position-reversed-nucleotide-content">
                    <div id="reversed-position-nucleotide-content-chart"></div>
                </div>

                <h3>Per Sequence GC Content</h3>
                <div class="chart-container" id="per-sequence-gc-content">
                    <div id="sequence-gc-content-chart"></div>
                </div>

            </section> <!-- END OF PER SEQUENCE DESCRIPTIVE STATS SECTION -->

        </div> <!-- END OF THE MAIN CONTENT AREA -->
    </div> <!-- END OF CONTAINER AREA -->

    
<!-- JAVASCRIPT TO POPULATE DATA AND CHARTS -->

    <script>

    </script>

</body>
</html>
"""

def put_file_details(html_template, filename):
    """
    Populates the placeholders {{filename}} and {{date}} in the HTML template.

    Args:
        html_template (str): The HTML template as a string.
        filename (str): The name of the file to insert into the template.

    Returns:
        str: The updated HTML template with placeholders replaced.
    """
    # Replace {{filename}} with the stripped filename
    html_template = html_template.replace("{{filename}}", filename)

    # Replace {{date}} with the current date and time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_template = html_template.replace("{{date}}", current_time)

    return html_template

def put_data(html_template, placeholder, data):
    """
    Replaces all occurrences of a placeholder in the HTML template with the provided data.

    Args:
        html_template (str): The HTML template as a string.
        placeholder (str): The placeholder to replace (e.g., "{{placeholder}}").
        data (str): The data to replace the placeholder with.

    Returns:
        str: The updated HTML template with placeholders replaced.

    Raises:
        ValueError: If the placeholder is not found in the HTML template.
    """
    if placeholder not in html_template:
        raise ValueError(f"Placeholder not found: {placeholder}")

    # Replace all occurrences of the placeholder with the data
    return html_template.replace(placeholder, data)

def escape_str(s):
    """
    Add \" around the string to escape it for HTML.
    """
    return '"' + s + '"'

def get_dataset_html_template(d1_stats, d2_stats, results):
    """
    Returns the HTML template for the report.
    """
    html_template = HTML_TEMPLATE

    # Replace placeholders with basic statistics for dataset 1
    html_template = put_data(html_template, "{{d1-filename}}", str(d1_stats['Filename']))
    html_template = put_data(html_template, "{{d1-number_of_sequences}}", str(d1_stats['Number of sequences']))
    html_template = put_data(html_template, "{{d1-number_of_bases}}", str(d1_stats['Number of bases']))
    html_template = put_data(html_template, "{{d1-unique_bases}}", ', '.join(x for x in d1_stats['Unique bases']))
    html_template = put_data(html_template, "{{d1-gc_content}}", f"{(d1_stats['%GC content']*100):.2f}%")
    html_template = put_data(html_template, "{{d1-dedup_sequences}}", str(d1_stats['Number of sequences left after deduplication']))


    # Replace placeholders with basic statistics for dataset 2
    html_template = put_data(html_template, "{{d2-filename}}", str(d2_stats['Filename']))
    html_template = put_data(html_template, "{{d2-number_of_sequences}}", str(d2_stats['Number of sequences']))
    html_template = put_data(html_template, "{{d2-number_of_bases}}", str(d2_stats['Number of bases']))
    html_template = put_data(html_template, "{{d2-unique_bases}}", ', '.join(x for x in d2_stats['Unique bases']))
    html_template = put_data(html_template, "{{d2-gc_content}}", f"{(d2_stats['%GC content']*100):.2f}%")
    html_template = put_data(html_template, "{{d2-dedup_sequences}}", str(d2_stats['Number of sequences left after deduplication']))    

    return html_template