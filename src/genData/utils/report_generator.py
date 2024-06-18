import os
import json
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def generate_html_report(stats_dict, results, output_path):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """
    pass


def generate_text_report(stats_dict, results, output_path):
    with open(output_path, 'w') as file:
        file.write('Sequence Statistics Report\n')
        for key, value in stats_dict.items():
            file.write(f'\n{key}: {results[key]}\n')
            if isinstance(value, dict):
                for sub_key, sub_value in value.items():
                    file.write(f'  {sub_key}: {sub_value}\n')
            else:
                file.write(f'  {value}\n')

def generate_simple_report(results, output_path):
    with open(output_path, 'w') as file:
        for key, value in results.items():
            if key == 'Filename':
                continue
            file.write(f'{key}: {value}\n')
    