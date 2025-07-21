HTML_CLUSTER_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Similar Sequences Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
        }
        h1 {
            text-align: center;
        }
        .cluster {
            border: 1px solid #ccc;
            margin-bottom: 20px;
            padding: 15px;
            border-radius: 5px;
        }
        .cluster h2 {
            margin-top: 0;
        }
        .section-title {
            font-weight: bold;
            margin-top: 10px;
        }
        pre {
            background: #f9f9f9;
            padding: 10px;
            overflow-x: auto;
        }
    </style>
</head>
<body>

<h1>Similar Sequences Found in Train vs Test Dataset</h1>

#clusters

</body>
</html>
"""

def get_train_test_html_template(clusters):
    if not clusters:
        return HTML_CLUSTER_TEMPLATE.replace("#clusters", "<h2>No similar sequences found.</h2>")

    cluster_blocks = []

    for cluster in clusters:
        cluster_html = f"""
        <div class="cluster">
            <h2>Cluster #{cluster['cluster']}</h2>
            <div class="section-title">Train Sequences:</div>
            <pre>{chr(10).join(cluster.get('train', []))}</pre>
            <div class="section-title">Test Sequences:</div>
            <pre>{chr(10).join(cluster.get('test', []))}</pre>
        </div>
        """
        cluster_blocks.append(cluster_html)

    return HTML_CLUSTER_TEMPLATE.replace("#clusters", "\n".join(cluster_blocks))