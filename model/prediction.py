import pandas as pd

def parse_contacts_to_dict(contact_file):
    '''
    For setting up restraint file from a getcontacts file.
    
    '''
    # TODO make it consider frequency and use that as confidence

    
    # Initialize the dictionary to store contact information
    contact_dict = {}

    with open(contact_file, 'r') as f:
        for line in f.readlines():
            if line.startswith("#") or not line.strip():
                # Skip comments or empty lines
                continue

            frame, ctype, resa, resb, dist = line.split()
            resa =  ":".join(resa.split(":")[:-1])
            resb =  ":".join(resb.split(":")[:-1])

            # Parse contact pair and distance
            contact_pair = tuple(sorted((resa, resb)))  # Ensure consistent ordering
            distance = float(dist)

            # Check if the contact pair already exists in the dictionary
            if contact_pair not in contact_dict:
                # Initialize with the current distance as both min and max
                contact_dict[contact_pair] = {
                    "min_distance": distance,
                    "max_distance": distance
                }
            else:
                # Update min and max distances as needed
                contact_dict[contact_pair]["min_distance"] = min(
                    contact_dict[contact_pair]["min_distance"], distance
                )
                contact_dict[contact_pair]["max_distance"] = max(
                    contact_dict[contact_pair]["max_distance"], distance
                )

    return contact_dict

def contact_dict_to_dataframe(contact_dict, resid_range=None, keep_chain=None):
    '''
    contact_dict : from parse_contacts_to_dict()

    resid_range : tuple
        minimum and maximum residue_id to include.
        e.g. only include residues 10 to 100 == (10, 100)

    keep_chain : str
        A chain ID if you only want to include results from a single chain
    
    Returns
    -------
    pd.DataFrame
    
    Can export to_csv(index=False) to use as a chai-1 restraint file.
    https://github.com/chaidiscovery/chai-lab/tree/main/examples/restraints
    '''
    restraints = []

    for idx, (contact_pair, distances) in enumerate(contact_dict.items()):
        # Extract chain and residue info
        resa, resb = contact_pair
        chA, resA, res_idxA = resa.split(':')[0], resa.split(':')[1], resa.split(':')[2]
        chB, resB, res_idxB = resb.split(':')[0], resb.split(':')[1], resb.split(':')[2]
        try:
            resA = convert_aa_code(resA)
        except:
            resA = "X"
        try:
            resB = convert_aa_code(resB)
        except:
            resB = "X"
            
        # Filter based on residue range
        if resid_range is not None:
            if not (resid_range[0] <= int(res_idxA) <= resid_range[1] and resid_range[0] <= int(res_idxB) <= resid_range[1]):
                continue

        # Filter based on chain
        if keep_chain is not None:
            if not (chA == keep_chain and chB == keep_chain):
                continue

       
        # Append the restraint row as a dictionary
        restraints.append({
            "restraint_id": f"restraint{idx}",
            "chainA": chA,
            "res_idxA": f"{resA}{res_idxA}",
            "chainB": chB,
            "res_idxB": f"{resB}{res_idxB}",
            "connection_type": "contact",
            "confidence": 1.0,
            "min_distance_angstrom": distances["min_distance"],
            "max_distance_angstrom": distances["max_distance"],
            "comment": f"Restraint for {resa}-{resb}"
        })

    # Convert the list of restraints to a pandas DataFrame
    df = pd.DataFrame(restraints)
    return df
