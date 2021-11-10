import sys
import json
import os
import networkx as nx
import matplotlib.pyplot as plt
import multiprocessing

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('branching_rules_input_file_path')
parser.add_argument('export_dir_path')
parser.add_argument('--export_first_n_branching_rules', type=int)

args = parser.parse_args()

if not os.path.exists(args.export_dir_path):
    os.makedirs(args.export_dir_path)

with open(args.branching_rules_input_file_path, 'r') as f:
    branching_rules = [json.loads(line) for line in f]
    for i, br in enumerate(branching_rules):
        br['file_idx'] = i+1

if args.export_first_n_branching_rules is not None:
    branching_rules = branching_rules[:args.export_first_n_branching_rules]

export_dir_path = args.export_dir_path
images_dir_path = export_dir_path + '/images'

for d in [export_dir_path, images_dir_path]:
    os.makedirs(d, exist_ok=True)


def generate_br_image(br):
    print('generating branching_rule={!r}'.format(br['file_idx']))
    plt.clf()

    graph = nx.Graph()
    graph.add_nodes_from(range(br['graph']['n']))

    node_color = []
    for v in range(br['graph']['n']):
        node_color.append('red' if v in br['red_vertices'] else 'lightblue')

    graph.add_edges_from(br['graph']['edges'])

    nx.draw_spring(graph, with_labels=True, node_size=1000, width=2, font_size=20, node_color=node_color)

    file_idx = br['file_idx']
    plt.savefig(f'{images_dir_path}/graph{file_idx}.png')

with multiprocessing.Pool(8) as p:
    p.map(generate_br_image, branching_rules)


def split(l, cnt):
    ll = []
    for e in l:
        ll.append(e)
        if len(ll) >= cnt:
            yield ll
            ll = []
    if ll:
        yield ll

generate_html_branches_cols = 3

def generate_html_branches(branches):
    branches = sorted([sorted(b) for b in branches])

    branches_len = (len(branches) + generate_html_branches_cols - 1) // generate_html_branches_cols * generate_html_branches_cols
    branches = (branches + [[] for _ in range(generate_html_branches_cols)])[:branches_len]

    html_output = ''
    for branches_row in split(branches, generate_html_branches_cols):
        html_output += '<tr>'
        for branch in branches_row:
            html_output += '<td>{}</td>'.format(', '.join([str(b) for b in branch]))
        html_output += '</tr>'

    return html_output

def generate_html_branching_rule(br):
    html_branching_rule_template = r'''
    <table class='branchingRuleImageTable'><tr>
    <td><img width="300px" src="{br_image_path}"/></td></tr><tr>
    <td><table width="300px" class='branchingRuleTable'>
        <tr><td colspan="3" class='branchingRuleNumber'># {br_index}</td></tr>
        <tr><td colspan="3" class='branchingRuleTableHeader'>Type:</td></tr>
        <tr><td colspan="3">{br_type}</td></tr>
        <tr><td colspan="3" class='branchingRuleTableHeader'>Branching factor:</td></tr>
        <tr><td colspan="3">{br_bf}</td></tr>
        <tr><td colspan="3" class='branchingRuleTableHeader'>Minimal branches:</td></tr>
        {br_minimal_branches}
        <tr><td colspan="3" class='branchingRuleTableHeader'>Dominance-free branches:</td></tr>
        {br_dominance_free_branches}
        <tr><td colspan="3" class='branchingRuleTableHeader'>Final branches:</td></tr>
        {br_final_branches}
    </table></td>
    </tr></table>
    '''.format(
        br_image_path = 'images/graph{}.png'.format(br['file_idx']),
        br_index = br['file_idx'],
        br_type = br['type'],
        br_bf = br['bf'],
        br_minimal_branches = generate_html_branches(br['minimal_branches']),
        br_dominance_free_branches = generate_html_branches(br['dominance_free_branches']),
        br_final_branches = generate_html_branches(br['branches']),
    )

    return html_branching_rule_template

generate_html_branching_rules_page_col = 6

def generate_html_branching_rules_page(branching_rules):
    html_output = ''
    html_output += '<table>\n'

    branching_rules_len = (len(branching_rules) + generate_html_branching_rules_page_col - 1) // generate_html_branching_rules_page_col * generate_html_branching_rules_page_col
    branching_rules = (branching_rules + [None for _ in range(generate_html_branching_rules_page_col)] )[:branching_rules_len]

    for branching_rules_row in split(branching_rules, generate_html_branching_rules_page_col):
        html_output += '<tr>\n'

        for branching_rule in branching_rules_row:
            html_output += '<td>\n'
            if branching_rule:
                html_output += generate_html_branching_rule(branching_rule)
            html_output += '</td>\n'

        html_output += '</tr>\n'

    html_output += '</table>'

    return html_output

def generate_html_index_page(branching_rules):
    html_style = r'''
<style>
    .branchingRuleImageTable, .branchingRuleImageTable td {
        border: 1px solid black;
        border-collapse: collapse;
        padding: 0px;
    }
    .branchingRuleTable, .branchingRuleTable td {
        border-collapse: collapse;
        border: 1px solid black;
        padding: 3px;
    }
    tr {
        vertical-align: top;
    }
    .branchingRuleTableHeader {
        background-color: orange;
    }
    .branchingRuleNumber {
        background-color: red;
    }

</style>
'''
    html_index_page_template = r'''
<!DOCTYPE html>
<html>
<head>
<title></title>
{style}
</head>
<body>
{html}
</body>
</html>
'''.format(style = html_style, html = generate_html_branching_rules_page(branching_rules))
    return html_index_page_template

with open(export_dir_path + '/index.html', 'w') as f:
    f.write(generate_html_index_page(branching_rules))

