import json

def convert_ints_to_strings(obj):
    """Recursively converts integer values in a dictionary or list to strings."""
    if isinstance(obj, int):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: convert_ints_to_strings(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_ints_to_strings(elem) for elem in obj]
    return obj

# Specify the input and output filenames
input_filename = "ecdlpcurvesnaive_extra.json"  # Save your JSON data to this file
output_filename = "ecdlpcurvesnaive_extra.json"

try:
    # Load the JSON data from the input file
    with open(input_filename, "r") as f:
        loaded_data = json.load(f)

    # Convert integer values to strings
    modified_data = convert_ints_to_strings(loaded_data)

    # Save the modified data to the output JSON file
    with open(output_filename, "w") as f:
        json.dump(modified_data, f, indent=4)

    print(f"Successfully loaded '{input_filename}', converted integers to strings, and saved to '{output_filename}'")

except FileNotFoundError:
    print(f"Error: Input file '{input_filename}' not found.")
except json.JSONDecodeError:
    print(f"Error: Could not decode JSON from '{input_filename}'. Please ensure it's valid JSON.")
except Exception as e:
    print(f"An error occurred: {e}")