import os
import sys
import pandas as pd
import plotly.graph_objects
from ._version import get_versions
from .get_scores import chim_score
from pkg_resources import resource_filename


def reshape_tax_levels(df, tax_levels):
    """Reshape tax data.

    Each col of df is a tax_level, each pair of adjacent levels is added to
    the return df in the form of:
         1) item in first tax level col (source)
         2) item in second tax level col (target)
         3) taxlevel of item in first col (tax level of source)

    Arguments:
        df (pandas.DataFrame): With each column being a taxlevel and rows being contigs.
        tax_levels (list): Of tax levels to pair together

    Returns:
        pandas.DataFrame: Reshaped to the form of source,target,source_tax_level
    """
    data = []
    for pair in zip(tax_levels, tax_levels[1:]):
        data += [row + [pair[0]] for row in df[list(pair)].values.tolist()]
    return pd.pandas.DataFrame(data, columns=['source', 'target', 'source_tax_level'])


def create_cat_codes_from_df(df):
    """Create category codes for each unique value in a df.

    Creates a dict mapping each unique value in a df to a number.

    Arguments:
        df (pandas.DataFrame): Containing values that need unique codes.

    Returns:
        Dict:Value: unique number
    """
    return {x: i for i, x in enumerate(set(df.values.ravel('K')))}


def convert_data(data, ref_dict):
    """Convert data using a feference dict.

    Replace all values in data using mapping in ref_dict.
    Data can be a pandas.DataFrame or iteable.

    Arguments:
        data (iter or pandas.DataFrame): Data to be converted
        ref_dict (Dict): Mapping of values to replace in data

    Returns:
        iter or pandas.DataFrame: data with values replaced
    """
    if isinstance(data, pd.pandas.DataFrame):
        return data.replace(ref_dict)
    else:
        return [ref_dict.get(item, item) for item in data]


def group_identical_rows(df):
    """Group identical rows.

    Merges all identical rows in a pandas.DataFrame and adds a col
    with the original count of rows.

    Arguments:
        df (pandas.DataFrame): With rows to be grouped.

    Returns:
        pandas.DataFrame: With duplicate rows merged and count col added.
    """
    columns = df.columns.tolist()
    return df.groupby(columns).size().to_frame('count').reset_index()


def extract_node_data(base_data, cat_codes):
    """Get node data for sankey plot.

    Data needs to be of the form:
        1) node: Int of each node in diagram
        2) colour: Hexcode of the node colour
        3) label: the label to assign to each node

    Arguments:
        base_data (pandas.DataFrame): With each column being a taxlevel and rows being contigs.
        cat_codes (Dict): Mapping of each clade in base_data to unique number.

    Returns:
        pandas.DataFrame: with node,colour,label columns
    """
    node_colours = {'kingdom': '#50514f',
                    'phylum': '#f25f5c',
                    'family': '#ffe066',
                    'genus': '#92AE83',
                    'species': '#78A1BB',
                    'contig': '#86BBBD'}
    nodes = list(cat_codes.keys())
    colour_dict = {}
    label_dict = {}
    for level in base_data.columns:
        for item in base_data[level].unique():
            colour_dict[item] = node_colours[level]
            if level != 'contig':
                label_dict[item] = item
    node_colours = [colour_dict.get(x, 'black') for x in nodes]
    node_labels = [label_dict.get(x, '') for x in nodes]
    return pd.pandas.DataFrame(list(zip(nodes, node_colours, node_labels)),
                               columns=['node', 'colour', 'label'])


def prepare_data(tax_data, tax_levels):
    """Prepare all data needed for sankey plot.

    Prepares node data and link data.

    Arguments:
        tax_data (pandas.DataFrame): With each column being a taxlevel and rows being contigs.
        tax_levels (list): Of tax levels to consider

    Returns:
        pandas.DataFrame: node_data
        pandas.DataFrame: link_data
    """

    node_colours = {'kingdom': '#50514f',
                    'phylum': '#f25f5c',
                    'family': '#ffe066',
                    'genus': '#92AE83',
                    'species': '#78A1BB',
                    'contig': '#86BBBD'}
    link_colours = {'kingdom': 'rgba(80,81,79,0.4)',
                    'phylum': 'rgba(242,95,92,0.4)',
                    'family': 'rgba(255,224,102,0.4)',
                    'genus': 'rgba(146,174,131,0.4)',
                    'species': 'rgba(120,161,187,0.4)',
                    'contig': 'rgba(134,187,189,0.4'}
    base_data = reshape_tax_levels(tax_data, tax_levels)
    cat_codes = create_cat_codes_from_df(base_data)
    link_data = group_identical_rows(base_data)
    link_data['sourceID'] = convert_data(link_data['source'], cat_codes)
    link_data['targetID'] = convert_data(link_data['target'], cat_codes)
    link_data['node_colours'] = convert_data(link_data['source_tax_level'],
                                             node_colours)
    link_data['link_colours'] = convert_data(link_data['source_tax_level'],
                                             link_colours)
    node_data = extract_node_data(tax_data[tax_levels], cat_codes)
    return node_data, link_data


def prepare_plot_data(node_data, link_data):
    """Create plotly Figure instance.

    Arguments:
        node_data (pandas.DataFrame): node_data
        link_data (pandas.DataFrame): link_data

    Returns:
        plotly.graph_objects.Figure: Sankey plot.
    """
    plot_data = {
            "data": [
                {
                    "type": "sankey",
                    "orientation": "h",
                    "arrangement": "freeform",
                    "node": {
                        "pad": 5,
                        "thickness": 8,
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
                        "customdata": link_data['target'],
                        "hovertemplate": '%{value:.g} genes from '
                                         '[%{source.label}] are assigned to '
                                         '[%{customdata}].<br />'
                                         '<extra></extra>'
                    }
                }],
            "layout": {
                "margin": {"t": 0, "l": 0, "r": 5},
                "font": {
                    "size": 10
                },

            }
        }
    return plotly.graph_objects.Figure(plot_data)


def get_html_template():
    """Read in HTML template.

    Returns:
        str: HTML Template
    """
    template_path = resource_filename(__name__, 'data/template.html')
    with open(template_path, 'r') as f:
        return f.read()


def create_html(plot_data, genome_name, display_info, levels_info):
    """Compile final HTML output.

    Put the plot ond ohter data in to HTML template

    Arguments:
        plot_data (plotly.graph_object.Figure): sankey plot
        genome_name (str): Name of genome
        display_info (str): Showing how many contigs were used to produce plot
        levels_info (str): Showing Tax levels being shown in plot

    Returns:
        str: Complete HTML
    """
    return get_html_template().format(plot=plot_data.to_html(),
                                      genome_name=genome_name,
                                      display_info=display_info,
                                      levels_info=levels_info,
                                      version=get_versions()['version'])


def parse_tax_levels_arg(tax_levels):
    """Parse contig level argument sting.

    Need to convert comma seperated input string to list to be used later.

    Arguments:
        tax_levels (str): commanseperated tax_levels

    Returns:
        list: tax levels to be used in plot
    """
    tax_levels = [x.strip() for x in tax_levels.split(',')]
    allowed = ['kingdom', 'phylum', 'family', 'genus', 'species', 'contig']
    if len(tax_levels) < 2:
        sys.exit('[Error] Need to provide at least 2 tax_levels.')
    for tax_level in tax_levels:
        if tax_level not in allowed:
            sys.exit(f'[Error] {tax_level} not known.'
                     f'Allowed: {",".join(allowed)}')
    return tax_levels


def create_viz_from_diamond_file(diamond_file, gene_count, tax_levels,
                                 contig_display_num, remove_minor_clade_level):
    """Create sankey plot.

    Uses diamond plot as input.

    Arguments:
        diamond_file (str): GUNC diamond output file path
        gene_count (int): Count of genes in original fasta
        tax_levels (str): Commaseperated taxlevels to consider in plot
        contig_display_num (int): Number of contigs to use for plot
        remove_minor_clade_level (str): Tax level at which to remove minor clades

    Returns:
        str: HTML to write to disk
    """
    if 'gtdb' in os.path.basename(diamond_file):
        db = 'gtdb_95'
    else:
        db = 'progenomes_2.1'
    tax_data, genome_name, cutoff = chim_score(diamond_file,
                                               gene_count,
                                               db=db,
                                               plot=True)
    total_contigs = len(tax_data['contig'].unique())
    if total_contigs > contig_display_num:
        print(f'[INFO] Subsampling data to display {contig_display_num} contigs.')
        tax_data = tax_data.groupby(remove_minor_clade_level).filter(
            lambda x: len(x) > cutoff)
        top_contigs = tax_data['contig'].value_counts().head(contig_display_num).index
        tax_data = tax_data[tax_data['contig'].isin(top_contigs)]
    tax_levels = parse_tax_levels_arg(tax_levels)
    node_data, link_data = prepare_data(tax_data, tax_levels)
    viz_data = prepare_plot_data(node_data, link_data)
    if contig_display_num > total_contigs:
        contig_display_num = total_contigs
    display_info = f'Displaying data from {contig_display_num}/{total_contigs} contigs.'
    levels_info = f'{" > ".join(tax_levels)}'
    return create_html(viz_data, genome_name, display_info, levels_info)
