import pandas as pd
from ..visualisation import *


def test_create_cat_codes_from_df():
    df = pd.DataFrame(['a', 'b', 'c'])
    cat_codes = create_cat_codes_from_df(df)
    assert set(cat_codes.values()) == {0, 1, 2}
    assert set(cat_codes.keys()) == {'a', 'b', 'c'}


def test_convert_data():
    data = [1, 2, 3]
    ref_dict = {1: 'a', 2: 'b', 3: 'c'}
    result = ['a', 'b', 'c']
    assert convert_data(data, ref_dict) == result
    df_data = pd.DataFrame(data)
    pd.testing.assert_frame_equal(convert_data(df_data, ref_dict),
                                  pd.DataFrame(result))


def test_group_identical_rows():
    df = pd.DataFrame([1, 2, 3, 1], columns=['num'])
    result = pd.DataFrame([[1, 2], [2, 1], [3, 1]], columns=['num', 'count'])
    grouped = group_identical_rows(df)
    pd.testing.assert_frame_equal(grouped, result)
