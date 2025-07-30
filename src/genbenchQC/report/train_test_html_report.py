HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Similar Sequences Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: flex;
            max-width: 100%;
        }
        .sidebar {
            width: 250px;
            background: #f4f4f4;
            padding: 20px;
            box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
            height: 100vh;
            position: fixed;
            overflow-y: auto;
            z-index: 1000;
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
            width: calc(100% - 350px);
            overflow-x: hidden;
        }
        section {
            margin-bottom: 50px;
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
        .data-item {
            font-size: 1.2em;
            margin-bottom: 10px;
        }
        .data-item span {
            font-weight: bold;
            font-size: 1em;
        }
        .cluster {
            border: 1px solid #ccc;
            margin-bottom: 20px;
            padding: 15px;
            border-radius: 5px;
        }
        pre {
            background: #f9f9f9;
            padding: 10px;
            overflow-x: auto;
        }
    </style>
</head>
<body>

    <div class="sidebar">
        <h2>Navigation</h2>
        <a href="#basic-descriptive-statistics">Basic Statistics</a>
        <a href="#clusters-section">Clusters</a>
    </div>

    <div class="content">
        <h1>Similar Sequences Found in Train vs Test Dataset</h1>

        <section id="basic-descriptive-statistics">
            <h2>Basic Descriptive Statistics</h2>
            <div class="data-item"><span>Train set filename:</span> {{train_filename}}</div>
            <div class="data-item"><span>Number of sequences in train set:</span> {{number_of_sequences_train}}</div>
            <div class="data-item"><span>Test set filename:</span> {{test_filename}}</div>
            <div class="data-item"><span>Number of sequences in test set:</span> {{number_of_sequences_test}}</div>
        </section>

        <section id="clusters-section">
            <h2>Clusters of Similar Sequences</h2>
            {{clusters}}
        </section>

    </div>

</body>
</html>
"""

def get_train_test_html_template(clusters, filename_train, sequences_train, filename_test, sequences_test):

    html_template = HTML_TEMPLATE

    html_template = html_template.replace("{{train_filename}}", str(filename_train))
    html_template = html_template.replace("{{test_filename}}", str(filename_test))
    html_template = html_template.replace("{{number_of_sequences_train}}", str(len(sequences_train)))
    html_template = html_template.replace("{{number_of_sequences_test}}", str(len(sequences_test)))

    if not clusters:
        return html_template.replace("{{clusters}}", "<h2>No similar sequences found.</h2>")

    cluster_blocks = []
    max_seq_display = 1000
    n_sequences = 0

    for cluster in clusters:
        train_sequences = cluster.get('train', [])
        test_sequences = cluster.get('test', [])
        if n_sequences + len(train_sequences) > max_seq_display:
            n_sequences += len(train_sequences)
            train_sequences = train_sequences[:2] + ["..."] if len(train_sequences) > 2 else train_sequences
            test_sequences = test_sequences[:2] + ["..."] if len(test_sequences) > 2 else test_sequences
        elif n_sequences + len(test_sequences) > max_seq_display:
            n_sequences += len(train_sequences) + len(test_sequences)
            test_sequences = test_sequences[:2] + ["..."] if len(test_sequences) > 2 else test_sequences
        else:
            n_sequences += len(train_sequences) + len(test_sequences)

        cluster_html = f"""
        <div class="cluster">
            <h2>Cluster #{cluster['cluster']}</h2>
            <div class="section-title">Train Sequences:</div>
            <pre>{chr(10).join(train_sequences)}</pre>
            <div class="section-title">Test Sequences:</div>
            <pre>{chr(10).join(test_sequences)}</pre>
        </div>
        """
        cluster_blocks.append(cluster_html)
        if n_sequences >= max_seq_display:
            cluster_blocks = [f"<b>Note:</b> There are too many clusters to display ({len(clusters)} clusters). Showing only part of the sequences. If you want to access all the clusters, toggle 'json' format and refer to the json report.  <br><br/>"] + cluster_blocks
            break

    return html_template.replace("{{clusters}}", "\n".join(cluster_blocks))