# Expanded tables

## Expanded

This style is ideal for very large tables that must scroll horizontally. It gives them more space and adds a drop shadow to make the scrolling area clearer.

:::{list-table}
---
header-rows: 1
class: expanded
---

* - Toolkit
  - residue_name
  - residue_number
  - chain_id
* - PDB file ATOM/HETATM columns
  - Columns 18-20 (`resName`)
  - Columns 23-26 (`resSeq`)
  - Columns 22 (`chainID`)
* - PDBx/MMCIF fields
  - `label_comp_id`
  - `label_seq_id`
  - `label_asym_id`
* - OpenFF getter (defined)
  - ```
    atom.metadata["residue_name"]
    ```
  - ```
    atom.metadata["residue_number"]
    ```
  - ```
    atom.metadata["chain_id"]
    ```
* - OpenFF getter (undefined)
  - ```
    "residue_name" not in atom.metadata
    ```
  - ```
    "residue_number" not in atom.metadata
    ```
  - ```
    "chain_id" not in atom.metadata
    ```
* - OpenFF setter (defined)
  - Defined: `atom.metadata["residue_name"] = X`
  - ```
    atom.metadata["residue_number"] = X
    ```
  - ```
    atom.metadata["chain_id"] = X
    ```
* - OpenFF setter (undefined)
  - ```
    del atom.metadata['residue_name']
    ```
  - ```
    del atom.metadata['residue_number']
    ```
  - ```
    del atom.metadata['chain_id']
    ```
:::


## Not expanded

This style doesn't make it obvious that the table can scroll.

:::{list-table}
---
header-rows: 1
---

* - Toolkit
  - residue_name
  - residue_number
  - chain_id
* - PDB file ATOM/HETATM columns
  - Columns 18-20 (`resName`)
  - Columns 23-26 (`resSeq`)
  - Columns 22 (`chainID`)
* - PDBx/MMCIF fields
  - `label_comp_id`
  - `label_seq_id`
  - `label_asym_id`
* - OpenFF getter (defined)
  - ```
    atom.metadata["residue_name"]
    ```
  - ```
    atom.metadata["residue_number"]
    ```
  - ```
    atom.metadata["chain_id"]
    ```
* - OpenFF getter (undefined)
  - ```
    "residue_name" not in atom.metadata
    ```
  - ```
    "residue_number" not in atom.metadata
    ```
  - ```
    "chain_id" not in atom.metadata
    ```
* - OpenFF setter (defined)
  - Defined: `atom.metadata["residue_name"] = X`
  - ```
    atom.metadata["residue_number"] = X
    ```
  - ```
    atom.metadata["chain_id"] = X
    ```
* - OpenFF setter (undefined)
  - ```
    del atom.metadata['residue_name']
    ```
  - ```
    del atom.metadata['residue_number']
    ```
  - ```
    del atom.metadata['chain_id']
    ```
:::