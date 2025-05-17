import json

def quote_large_integers(data):
    if isinstance(data, dict):
        new_dict = {}
        for key, value in data.items():
            new_dict[key] = quote_large_integers(value)
        return new_dict
    elif isinstance(data, list):
        new_list = []
        for item in data:
            new_list.append(quote_large_integers(item))
        return new_list
    elif isinstance(data, int):
        # Convert integer to string and quote it
        return str(data)
    else:
        return data

def process_json_file(input_filename, output_filename):
    try:
        with open(input_filename, 'r') as infile:
            json_data = json.load(infile)

        modified_data = quote_large_integers(json_data)

        with open(output_filename, 'w') as outfile:
            json.dump(modified_data, outfile, indent=4)  # Use indent for better readability

        print(f"Successfully processed '{input_filename}' and saved the quoted integers to '{output_filename}'")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{input_filename}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    input_file = "shamir.json"  # Replace with the actual name of your input file
    output_file = "shamirnew.json" # Replace with the desired name for the output file
    process_json_file(input_file, output_file)