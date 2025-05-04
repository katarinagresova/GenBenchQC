from datetime import datetime

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HTML Report Output</title>
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
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

        section {
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
    <div class="container">
        <div class="sidebar">
            <h2>Navigation</h2>
            <a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a>
            <div style="margin-left: 15px;">
                <a href="#filename">Filename</a>
                <a href="#num-sequences">Number of sequences</a>
                <a href="#num-bases">Number of bases</a>
                <a href="#unique-bases">Unique bases</a>
                <a href="#gc-content">%GC content</a>
                <a href="#dedup-sequences">Number of sequences left after deduplication</a>
            </div>
            <a href="#general-descriptive-statistics">General Descriptive Statistics</a>
            <div style="margin-left: 15px;">
                <a href="#sequence-lengths">Sequence lengths</a>
                <a href="#sequence-duplication-levels">Sequence duplication levels</a>
            </div>
            <a href="#per-sequence-descriptive-stats">Per Sequence Descriptive Stats</a>
            <div style="margin-left: 15px;">
                <a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a>
                <a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a>
                <a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a>
                <a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a>
                <a href="#per-sequence-gc-content">Per Sequence GC Content</a>
            </div>
        </div>

        <div class="content">
        <h1>HTML Report Output</h1>

            <section id="basic-descriptive-statistics">
                <h2>Basic Descriptive Statistics</h2>
                <div class="data-item" id="filename">
                    <span>Filename:</span> <!-- Filename will be displayed here -->
                </div>
                <div class="data-item" id="num-sequences">
                    <span>Number of sequences:</span> <!-- Number of sequences will be displayed here -->
                </div>
                <div class="data-item" id="num-bases">
                    <span>Number of bases:</span> <!-- Number of bases will be displayed here -->
                </div>
                <div class="data-item" id="unique-bases">
                    <span>Unique bases:</span> <!-- Unique bases will be displayed here -->
                </div>
                <div class="data-item" id="gc-content">
                    <span>%GC content:</span> <!-- %GC content will be displayed here -->
                </div>
                <div class="data-item" id="dedup-sequences">
                    <span>Number of sequences left after deduplication:</span> <!-- Number of sequences left after deduplication will be displayed here -->
                </div>
            </section>

            <section id="general-descriptive-statistics">
                <h2>General Descriptive Statistics</h2>

                <h3>Sequence lengths</h3>
                <div class="chart-container" id="sequence-lengths">

                </div>
                <div id="sequence-duplication-levels">
                    <h3>Sequence duplication levels</h3>
                    <table>
                        <thead>
                            <tr>
                                <th class="count_column">Count</th>
                                <th class="sequence_column">Sequence</th>
                            </tr>
                        </thead>
                        <tbody>
                            <!-- Table rows will be dynamically populated -->
                        </tbody>
                    </table>
                </div>
            </section>

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

            </section>
        </div>
    </div>

    <script>
<!--Basic stats section-->
        var basicStats = {
            filename: {{filename}},
            numSequences: {{number_of_sequences}},
            numBases: {{number_of_bases}},
            uniqueBases: {{unique_bases}},
            gcContent: {{gc_content}},
            dedupSequences: {{dedup_sequences}}
        };

        // Populate the Basic Descriptive Statistics section with the example data
        document.getElementById("filename").innerHTML += basicStats.filename;
        document.getElementById("num-sequences").innerHTML += basicStats.numSequences;
        document.getElementById("num-bases").innerHTML += basicStats.numBases;
        document.getElementById("unique-bases").innerHTML += basicStats.uniqueBases.join(", ");
        document.getElementById("gc-content").innerHTML += basicStats.gcContent + "%";
        document.getElementById("dedup-sequences").innerHTML += basicStats.dedupSequences;

        var sequenceDuplicationLevels = {{sequence_duplication_levels}};

        var lengths = {{sequence_lengths}};

        // Building the Frequency Histogram
        var trace = {
            x: lengths,
            type: 'histogram',
        };
        var data = [trace];
        var layout = {
            xaxis: { title: 'Sequence Length' },
            yaxis: { title: 'Frequency' }
        };

        Plotly.newPlot('sequence-lengths', data, layout);


        // Populate table for sequence duplication levels
        var tableBody = document.querySelector("#sequence-duplication-levels tbody");
        for (var sequence in sequenceDuplicationLevels) {
            var row = document.createElement("tr");
            var countCell = document.createElement("td");
            var sequenceCell = document.createElement("td");

            countCell.textContent = sequenceDuplicationLevels[sequence];
            countCell.className = "count_column";
            sequenceCell.textContent = sequence;
            sequenceCell.className = "sequence_column";

            row.appendChild(countCell);
            row.appendChild(sequenceCell);
            tableBody.appendChild(row);
        }

        ///PER SEQUENCE PLOTS


        var perSequenceGCContent = {
            "seq1": 55.5,
            "seq2": 44.5
        };

        // Summary statistics for per sequence nucleotide content
        var nucleotideContentSummary = {{per_sequence_nucleotide_content_summary}};

        // Prepare data for the box plot
        var nucleotideData = [];
        for (var nucleotide in nucleotideContentSummary) {
            var stats = nucleotideContentSummary[nucleotide];
            nucleotideData.push({
                y: [stats.min, stats.q1, stats.median, stats.q3, stats.max],
                type: 'box',
                name: nucleotide,
                boxpoints: false // Disable individual data points
            });
        }

        // Layout for the plot
        var layoutNuc = {
            xaxis: { title: 'Nucleotide' },
            yaxis: { title: 'Counts' },
            showlegend: false
        };

        // Render the plot
        Plotly.newPlot('per-sequence-nucleotide-content', nucleotideData, layoutNuc);

        // Summary statistics for per sequence dinucleotide content
        var dinucleotideContentSummary = {{per_sequence_dinucleotide_content_summary}};

        // Prepare data for the box plot
        var dinucleotideData = [];
        for (var dinucleotide in dinucleotideContentSummary) {
            var stats = dinucleotideContentSummary[dinucleotide];
            dinucleotideData.push({
                y: [stats.min, stats.q1, stats.median, stats.q3, stats.max],
                type: 'box',
                name: dinucleotide,
                boxpoints: false // Disable individual data points
            });
        }

        // Layout for the plot
        var layoutNuc = {
            xaxis: { title: 'Dinucleotide' },
            yaxis: { title: 'Frequency' },
            showlegend: false
        };

        // Render the plot
        Plotly.newPlot('per-sequence-dinucleotide-content', dinucleotideData, layoutNuc);

        // Per position nucleotide content, plotted as lineplots
        var perPositionNucleotideContent = {{per_position_nucleotide_content}};
        var positions = Object.keys(perPositionNucleotideContent);
        var nucleotides = basicStats.uniqueBases;

        var scatterData = nucleotides.map(nucleotide => {
            return {
                y: positions.map(pos => perPositionNucleotideContent[pos][nucleotide] || 0),
                type: 'scatter',
                name: nucleotide,
            };
        });
        var scatterLayout = {
            xaxis: { title: 'Position' },
            yaxis: { title: 'Frequency' },
            showlegend: true,
        };
        Plotly.newPlot('per-position-nucleotide-content', scatterData, scatterLayout);

        // Reverse per position nucleotide content, plotted as lineplots
        var revPerPositionNucleotideContent = {{per_position_reversed_nucleotide_content}};
        var positions = Object.keys(revPerPositionNucleotideContent);
        var nucleotides = basicStats.uniqueBases;

        var scatterData = nucleotides.map(nucleotide => {
            return {
                y: positions.map(pos => revPerPositionNucleotideContent[pos][nucleotide] || 0),
                type: 'scatter',
                name: nucleotide,
            };
        });
        var scatterLayout = {
            xaxis: { title: 'Position' },
            yaxis: { title: 'Frequency' },
            showlegend: true,
        };
        Plotly.newPlot('per-position-reversed-nucleotide-content', scatterData, scatterLayout);

        var perSequenceGCContent = {{per_sequence_gc_content}};
        var gcContentData = {
            x: Object.values(perSequenceGCContent),
            type: 'histogram',
            name: 'GC Content'
        };
        var gcContentLayout = {
            xaxis: { title: 'GC Content (%)' },
            yaxis: { title: 'Counts' },
        };
        Plotly.newPlot('sequence-gc-content-chart', [gcContentData], gcContentLayout);

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

def get_html_template(stats):
    """
    Returns the HTML template for the report.
    """
    html_template = HTML_TEMPLATE

    # Replace placeholders with computed statistics
    html_template = put_data(html_template, "{{filename}}", escape_str(str(stats['Filename'])))
    html_template = put_data(html_template, "{{number_of_sequences}}", str(stats['Number of sequences']))
    html_template = put_data(html_template, "{{number_of_bases}}", str(stats['Number of bases']))
    html_template = put_data(html_template, "{{unique_bases}}", '[' + ', '.join(escape_str(x) for x in stats['Unique bases']) + ']')
    html_template = put_data(html_template, "{{gc_content}}", f"{(stats['%GC content']*100):.2f}")
    html_template = put_data(html_template, "{{dedup_sequences}}", str(stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{sequence_lengths}}", '[' + ', '.join(map(str, stats['Sequence lengths'].values())) + ']')
    html_template = put_data(html_template, "{{sequence_duplication_levels}}", str(stats['Sequence duplication levels']).replace("'", '"').replace(", ", ",\n"))
    html_template = put_data(html_template, "{{per_sequence_nucleotide_content_summary}}", str(stats['Nucleotide content summary']).replace("'", '"').replace(", ", ",\n"))
    html_template = put_data(html_template, "{{per_sequence_dinucleotide_content_summary}}", str(stats['Dinucleotide content summary']).replace("'", '"').replace(", ", ",\n"))
    # Take only first 100 positions from per position nucleotide content
    position_content = {pos: stats['Per position nucleotide content'][pos] for pos in list(stats['Per position nucleotide content'].keys())[:100]}
    html_template = put_data(html_template, "{{per_position_nucleotide_content}}", str(position_content).replace("'", '"').replace(", ", ",\n"))

    reverse_position_content = {pos: stats['Per position reversed nucleotide content'][pos] for pos in list(stats['Per position reversed nucleotide content'].keys())[:100]}
    html_template = put_data(html_template, "{{per_position_reversed_nucleotide_content}}", str(reverse_position_content).replace("'", '"').replace(", ", ",\n"))

    html_template = put_data(html_template, "{{per_sequence_gc_content}}", str(stats['Per sequence GC content']).replace("'", '"').replace(", ", ",\n"))

    return html_template