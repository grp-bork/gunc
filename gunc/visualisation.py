import plotly
import pandas as pd
from ._version import get_versions
from pkg_resources import resource_filename


def reshape_tax_levels(df, tax_levels):
    data = []
    for pair in zip(tax_levels, tax_levels[1:]):
        data += [row + [pair[0]] for row in df[list(pair)].values.tolist()]
    return pd.DataFrame(data, columns=['source', 'target', 'source_tax_level'])


def create_cat_codes_from_df(df):
    return {x: i for i, x in enumerate(set(df.values.ravel('K')))}


def convert_data(data, ref_dict):
    if isinstance(data, pd.DataFrame):
        return data.replace(ref_dict)
    else:
        return [ref_dict.get(item, item) for item in data]


def group_identical_rows(df):
    columns = df.columns.tolist()
    return df.groupby(columns).size().to_frame('count').reset_index()


def extract_node_data(base_data, cat_codes):
    nodes = list(cat_codes.keys())
    colours = ["#50514f", "#f25f5c", "#ffe066",
               "#92AE83", "#78A1BB", "#86BBBD"]
    colour_dict = {}
    label_dict = {}
    for level, colour in zip(base_data.columns, colours):
        for item in base_data[level].unique():
            colour_dict[item] = colour
            if level != 'contig':
                label_dict[item] = item
    node_colours = [colour_dict.get(x, 'black') for x in nodes]
    node_labels = [label_dict.get(x, '') for x in nodes]
    return pd.DataFrame(list(zip(nodes, node_colours, node_labels)),
                        columns=['node', 'colour', 'label'])


def prepare_data(tax_data):
    node_colours = {'kingdom': '#50514f',
                    'phylum': '#f25f5c',
                    'family': '#ffe066',
                    'genus': '#92AE83',
                    'specI': '#78A1BB',
                    'contig': '#86BBBD'}
    link_colours = {'kingdom': 'rgba(80,81,79,0.4)',
                    'phylum': 'rgba(242,95,92,0.4)',
                    'family': 'rgba(255,224,102,0.4)',
                    'genus': 'rgba(146,174,131,0.4)',
                    'specI': 'rgba(120,161,187,0.4)',
                    'contig': 'rgba(134,187,189,0.4'}
    tax_levels = ['kingdom', 'phylum', 'family', 'genus', 'specI', 'contig']
    base_data = reshape_tax_levels(tax_data, tax_levels)
    cat_codes = create_cat_codes_from_df(base_data)
    link_data = group_identical_rows(base_data)
    link_data['sourceID'] = convert_data(link_data['source'], cat_codes)
    link_data['targetID'] = convert_data(link_data['target'], cat_codes)
    link_data['node_colours'] = convert_data(link_data['source_tax_level'],
                                             node_colours)
    link_data['link_colours'] = convert_data(link_data['source_tax_level'],
                                             link_colours)
    node_data = extract_node_data(tax_data, cat_codes)
    return node_data, link_data


def prepare_plot_data(node_data, link_data, genome_name):
    plot_data = {
            "data": [
                {
                    "type": "sankey",
                    "orientation": "h",
                    "node": {
                        "pad": 5,
                        "thickness": 28,
                        "line": {
                            "color": "grey",
                            "width": 0.1
                        },
                        "label": node_data['label'],
                        "color": node_data['colour'],
                        "customdata": node_data['node'],
                        "hovertemplate": '[%{customdata}] has %{value:.g} '
                                         'genes assigned.<extra></extra>'
                    },
                    "link": {
                        "source": link_data['sourceID'],
                        "target": link_data['targetID'],
                        "value": link_data['count'],
                        "color": link_data['link_colours'],
                        "hovertemplate": '%{value:.g} genes from '
                                         '[%{source.label}] are assigned to '
                                         '[%{target.label}].<br />'
                                         '<extra></extra>'
                    }
                }],
            "layout": {

                "font": {
                    "size": 10
                },
                "updatemenus": [
                    {
                        "y": 1.1,
                        "x": 0.0,
                        "buttons": [
                            {
                                "label": "Light Theme  ",
                                "method": "relayout",
                                "args": ["paper_bgcolor", "white"]
                            },
                            {
                                "label": "Dark Theme",
                                "method": "relayout",
                                "args": ["paper_bgcolor", "black"]
                            }
                        ]
                    },
                    {
                        "y": 1.1,
                        "x": 0.1,
                        "buttons": [
                            {
                                "label": "Thick Nodes ",
                                "method": "restyle",
                                "args": ["node.thickness", 15]
                            },
                            {
                                "label": "Thin Nodes",
                                "method": "restyle",
                                "args": ["node.thickness", 8]
                            }
                        ]
                    },
                    {
                        "y": 1.1,
                        "x": 0.2,
                        "buttons": [
                            {
                                "label": "Fixed Move   ",
                                "method": "restyle",
                                "args": ["arrangement", "perpendicular"]
                            },
                            {
                                "label": "Free Move",
                                "method": "restyle",
                                "args": ["arrangement", "freeform"]
                            }
                        ]
                    },
                    {
                        "y": 1.1,
                        "x": 0.3,
                        "buttons": [
                            {
                                "label": "Horizontal   ",
                                "method": "restyle",
                                "args": ["orientation", "h"]
                            },
                            {
                                "label": "Vertical",
                                "method": "restyle",
                                "args": ["orientation", "v"]
                            }
                        ]
                    }
                ]
            }
        }
    return plotly.graph_objects.Figure(plot_data)


def create_data(tax_data, genome_name):
    node_data, link_data = prepare_data(tax_data['base_data'])
    plot_data = prepare_plot_data(node_data, link_data, genome_name)
    return plotly.io.to_json(plot_data, pretty=True)


def get_html_template():
    template_path = resource_filename(__name__, 'data/template.html')
    with open(template_path, 'r') as f:
        return f.read()


def create_html(plot_data, genome_name):
    plot = plotly.io.from_json(plot_data)
    return get_html_template().format(plot=plotly.graph_objects.Figure(plot),
                                      genome_name=genome_name,
                                      version=get_versions()['version'])


def write_html(plot, genome_name, out_file):
    out_html = get_html_template().format(plot=plot,
                                          genome_name=genome_name,
                                          version=get_versions()['version'])
    with open(out_file, 'w') as f:
        f.write(out_html)
