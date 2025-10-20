import re
import sys
import os
from collections import defaultdict

def extract_bioinformatics_details(log_file_path):
    """
    Parses a bioinformatics log file to extract tool versions and 
    associated database details from known patterns.
    """
    results = defaultdict(dict)
    
    # Regex to capture tool names and versions. It's robust for formats like:
    # 1. 'Tracked executable version: kneaddata v0.12.0'
    # 2. 'Tracked executable version: MetaPhlAn version 4.0.6 (1 Mar 2023)'
    TOOL_VERSION_REGEX = re.compile(
        r"Tracked executable version:\s*(?P<tool>[A-Za-z0-9]+.*?)[\s:]+(?:v|version\s+)?(?P<version>[\d\.\-]+(?:[\s\(][\s\w\d\.\-]+\))?)"
    )

    # Regex to capture database configuration variables common in biobakery workflows.
    CONFIG_DB_REGEX = re.compile(
        r"INFO: (?P<config_key>contaminate_databases|metaphlan_db|humann_databases|chocophlan_db|utility_mapping_db)\s*=\s*(?P<db_path>.*)"
    )
    
    try:
        with open(log_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # --- 1. Extract Tool Versions ---
                m_tool = TOOL_VERSION_REGEX.search(line)
                if m_tool:
                    tool_name = m_tool.group("tool").strip()
                    version = m_tool.group("version").strip()
                    
                    # Standardize tool name (case-insensitive check)
                    if "kneaddata" in tool_name.lower():
                        tool_name = "kneaddata"
                    elif "metaphlan" in tool_name.lower():
                        tool_name = "MetaPhlAn"
                    elif "humann" in tool_name.lower():
                        tool_name = "HUMAnN"
                    
                    # Store only the first version found for a tool
                    if not results[tool_name].get("Tool_Version"):
                         results[tool_name]["Tool_Version"] = version

                # --- 2. Extract Database Details from Config ---
                m_db_config = CONFIG_DB_REGEX.search(line)
                if m_db_config:
                    config_key = m_db_config.group("config_key")
                    db_path = m_db_config.group("db_path").strip()
                    
                    # Map config key to the corresponding tool
                    tool_map = {
                        "contaminate_databases": "kneaddata",
                        "metaphlan_db": "MetaPhlAn",
                        "humann_databases": "HUMAnN",
                        "chocophlan_db": "HUMAnN",
                        "utility_mapping_db": "HUMAnN",
                    }
                    
                    tool_name = tool_map.get(config_key, "Other_Databases")

                    # Use the last folder name as a guess for the database name
                    db_name_guess = os.path.basename(os.path.normpath(db_path))
                    
                    if db_path and db_name_guess:
                        db_entry = {
                            "Name_from_Path": db_name_guess,
                            "Full_Path": db_path,
                            "Config_Key": config_key
                        }
                        
                        # Initialize Database_Details structure
                        if "Database_Details" not in results[tool_name]:
                             results[tool_name]["Database_Details"] = {}
                        
                        # Use the config key as the database identifier to allow multiple DBs per tool
                        if config_key not in results[tool_name]["Database_Details"]:
                            results[tool_name]["Database_Details"][config_key] = db_entry

    except FileNotFoundError:
        print(f"Error: The file '{log_file_path}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        sys.exit(1)
        
    return results

def format_output(results):
    """Formats the extracted details into a readable string."""
    output = ["\n--- Extracted Bioinformatics Tool and Database Details ---"]
    
    # Define preferred order
    tools_order = ["kneaddata", "MetaPhlAn", "HUMAnN"]
    
    # Sort keys for consistent output, placing requested tools first
    sorted_tools = tools_order + [k for k in results if k not in tools_order]
    
    for tool_name in sorted_tools:
        if tool_name not in results:
            continue
            
        details = results[tool_name]
        version = details.get("Tool_Version", "Not Found")
        
        output.append(f"\n[ Tool: {tool_name} ]")
        output.append(f"  Version: {version}")
        
        db_details = details.get("Database_Details")
        if db_details:
            output.append("  Databases:")
            for db_key, db_info in db_details.items():
                db_key_clean = db_key.replace('_', ' ').title()
                output.append(f"    - Configuration Key: {db_key_clean}")
                output.append(f"      Name (Guessed from path): {db_info['Name_from_Path']}")
                output.append(f"      Full Path: {db_info['Full_Path']}")
        else:
            output.append("  Databases: No specific database details extracted.")
            
    return "\n".join(output)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_tool_details.py <path_to_log_file>")
        print("\nExample: python extract_tool_details.py anadama.log")
        sys.exit(1)

    log_file_path = sys.argv[1]
    
    extracted_data = extract_bioinformatics_details(log_file_path)
    
    # Print the final formatted results
    print(format_output(extracted_data))


    --- Extracted Bioinformatics Tool and Database Details ---

# Example Output
# [ Tool: kneaddata ]
#   Version: 0.12.0
#   Databases:
#     - Configuration Key: Contaminate Databases
#       Name (Guessed from path): hg37_and_human_contamination
#       Full Path: /n/huttenhower_lab/data/kneaddata_databases/hg37_and_human_contamination/

# [ Tool: MetaPhlAn ]
#   Version: 4.0.6
#   Databases: No specific database details extracted.

# [ Tool: HUMAnN ]
#   Version: 3.7
#   Databases: No specific database details extracted.