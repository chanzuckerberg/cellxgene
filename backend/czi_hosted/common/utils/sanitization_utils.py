import re


def sanitize_values_in_list(list_of_keys: list):
    """
    Returns a dictionary mapping of the old keys in the list of `list_of_keys` to its new, clean name that is both
    safe and unique.
    """

    if not all([isinstance(key, str) for key in list_of_keys]):
        raise Exception("List of keys to sanitize must contain all strings.")

    # Mask out [~/.] and anything outside the ASCII range.
    mask = re.compile(r"[^ -\-0-\[\]-\}]")
    clean_keys_list = [mask.sub("_", key) for key in list_of_keys]

    # Dedupe the clean keys list
    deduped_clean_keys_list = []
    for index, clean_key in enumerate(clean_keys_list):
        total_occurrences_of_clean_key = clean_keys_list.count(clean_key)
        total_occurrences_up_until_current_index = clean_keys_list[:index].count(clean_key)
        deduped_clean_keys_list.append(
            clean_key + "_" + str(total_occurrences_up_until_current_index + 1)
            if total_occurrences_of_clean_key > 1
            else clean_key
        )

    return dict(zip(list_of_keys, deduped_clean_keys_list))


def sanitize_keys_in_dictionary(dict_to_sanitize: dict):
    """
    Clean and dedupe the keys in the given dictionary.
    """

    clean_keys = sanitize_values_in_list(dict_to_sanitize.keys())
    for original_key, sanitized_key in clean_keys.items():
        if original_key != sanitized_key:
            dict_to_sanitize[sanitized_key] = dict_to_sanitize[original_key]
            del dict_to_sanitize[original_key]
