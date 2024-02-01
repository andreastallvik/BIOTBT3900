"""Functions that can be useful for dealing with the exp. data from the case study."""


def remove_neg_vals(df, inplace=False):
    """Remove negative values from a dataframe, for columns containing float."""    

    return_df = df.applymap(lambda x: 0 if isinstance(x, float) and x < 0 else x)

    if inplace:
        df = return_df
    else:
        return return_df