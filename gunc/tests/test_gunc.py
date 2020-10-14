import pytest
import pandas as pd
from ..get_scores import *
from pkg_resources import resource_filename

diamond_output = resource_filename(__name__,
                                   'test_data/tiny_test.diamond.out')
diamond_df = read_diamond_output(diamond_output)
ref_base_data_path = resource_filename(__name__,
                                       'test_data/tiny_test.base_data')
ref_base_data = pd.read_csv(ref_base_data_path)
