from typing import List


def split_list(input_list: List, num_lists: int) -> List[List]:
    """
    From a source list, creates a list of lists of the same elements of the original list,
    where each inner list is the length of the specified number. If the length of the original
    list is not evenly divided by the specified number, the last inner list will have a length
    less than this number.

    This provides a similar function as `itertools.batched()`, which
    was added in Python 3.12.

    Args:
        input_list: List of objects.
        num_lists: Number of lists to create.

    Returns: List of lists.

    """
    # find the desired number elements in each list
    list_len = len(input_list) // num_lists

    # create the first num_lists-1 lists
    nested_list = [
        input_list[index * list_len: (index + 1) * list_len] for index in range(num_lists - 1)
    ]
    # create the final list, including the division remainder number of elements
    nested_list.append(input_list[(num_lists - 1) * list_len:])
    return nested_list
