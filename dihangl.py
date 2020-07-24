import os

def master_mmCIF_renumber_function(input_mmCIF_files_were_found, default_input_path_to_mmCIF,
                                   default_input_path_to_SIFTS, default_output_path_to_mmCIF,
                                   default_mmCIF_num, gzip_mode):
    stepik = 0
    if not os.path.exists(default_output_path_to_mmCIF):
        os.makedirs(default_output_path_to_mmCIF)

    for mmCIF_ID in tqdm(input_mmCIF_files_were_found, leave=True, position=0, unit="files", desc="Renumbering mmCIF"):
        # print("Start Processing:", mmCIF_ID[:4])
        SIFTS_ID = mmCIF_ID[:4] + ".xml.gz"
        stepik = stepik + 1

        try:
            handle_SIFTS = gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_ID), 'rt')
        except FileNotFoundError:
            # print("No corresponding SIFTS for:", mmCIF_ID)

            mmcif_dict = 0
            for _ in range(3):
                try:
                    mmcif_dict = MMCIF2Dict.MMCIF2Dict(gzip.open(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_ID), 'rt'))
                    break
                except EOFError:
                    os.remove(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_ID))
                    download_master("mmCIF", [mmCIF_ID])
                    # print("Cannot open this file:", mmCIF_ID)
            if mmcif_dict == 0:
                return None

            os.chdir(default_output_path_to_mmCIF)
            io = MMCIFIO()
            io.set_dict(mmcif_dict)
            io.save(mmCIF_ID[:4] + "_no_SIFTS_out.cif")
            if gzip_mode == 0:
                gzip_for_output_files(mmCIF_ID[:4] + "_no_SIFTS_out.cif", gzip_mode)
                os.remove(mmCIF_ID[:4] + "_no_SIFTS_out.cif")
            os.chdir(currentDirectory)
            continue

        product_tree_SIFTS = 0
        for _ in range(3):
            try:
                product_tree_SIFTS = SIFTS_tree_parser(handle_SIFTS)
                break
            except EOFError:
                os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_ID))
                download_master("SIFTS", [mmCIF_ID])
                # print("Cannot open this file:", SIFTS_ID)
        if product_tree_SIFTS == 0:
            return None

        tuple_PDBe_for_PDB_and_tuple_PDB = product_tree_SIFTS[0]
        tuple_PDBe_for_UniProt_and_tuple_UniProt = product_tree_SIFTS[1]

        if tuple_PDBe_for_UniProt_and_tuple_UniProt == list():
            mmcif_dict = 0
            for _ in range(3):
                try:
                    mmcif_dict = MMCIF2Dict.MMCIF2Dict(
                        gzip.open(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_ID), 'rt'))
                    break
                except EOFError:
                    os.remove(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_ID))
                    download_master("mmCIF", [mmCIF_ID])
                    # print("Cannot open this file:", mmCIF_ID)
            if mmcif_dict == 0:
                return None
            # print("No UniProt in SIFTS!", mmCIF_ID[:4])
            os.chdir(default_output_path_to_mmCIF)
            io = MMCIFIO()
            io.set_dict(mmcif_dict)
            io.save(mmCIF_ID[:4] + "_no_UniProt_in_SIFTS_out.cif")
            if gzip_mode == 0:
                gzip_for_output_files(mmCIF_ID[:4] + "_no_UniProt_in_SIFTS_out.cif", gzip_mode)
                os.remove(mmCIF_ID[:4] + "_no_UniProt_in_SIFTS_out.cif")
            os.chdir(currentDirectory)

        product_of_SIFTS_data_parser = SIFTS_data_parser_for_mmCIF(tuple_PDBe_for_PDB_and_tuple_PDB,
                                                                   tuple_PDBe_for_UniProt_and_tuple_UniProt,
                                                                   default_mmCIF_num)
        df_PDBe_PDB_UniProt_without_null_index_PDBe = product_of_SIFTS_data_parser[0]
        df_PDBe_PDB_UniProt = product_of_SIFTS_data_parser[1]

        product_of_mmCIF_parser = mmCIF_parser(mmCIF_ID, default_input_path_to_mmCIF,
                                               df_PDBe_PDB_UniProt_without_null_index_PDBe, default_mmCIF_num)
        df_final_dropped_dup = product_of_mmCIF_parser[0]
        mmcif_dict = product_of_mmCIF_parser[1]
        _pdbx_poly_seq_scheme_auth_seq_num_before_change = product_of_mmCIF_parser[2]
        _atom_site_label_comp_id_list = product_of_mmCIF_parser[3]

        mmcif_dict_keys = mmcif_dict.keys()
        formed_columns = column_formation(mmcif_dict_keys)
        for n in tqdm(formed_columns, leave=False, desc="Renumbering " + mmCIF_ID[:4]):
            PDB_ins_code = n[3]
            auth_asym_id = n[2]
            auth_comp_id = n[1]
            auth_seq_id = n[0]
            renumber_small_tables(auth_comp_id, auth_asym_id, auth_seq_id, PDB_ins_code, mmcif_dict,
                                  df_final_dropped_dup, default_mmCIF_num)

        renum_struct_site_details(mmcif_dict, _atom_site_label_comp_id_list, df_final_dropped_dup, default_mmCIF_num)
        renum_pdbx_unobs_or_zero_occ_residues_auth_seq_id(mmcif_dict, df_PDBe_PDB_UniProt, default_mmCIF_num)
        renum_pdbx_poly_seq_scheme_auth_seq_num(mmcif_dict, df_final_dropped_dup, default_mmCIF_num)
        renum_pdbx_nonpoly_scheme_auth_seq_num(mmcif_dict, df_final_dropped_dup, default_mmCIF_num)

        # try:
        # _pdbx_poly_seq_scheme_auth_seq_num_after_change = mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"]
        # same_len = len([i for i, j in zip(_pdbx_poly_seq_scheme_auth_seq_num_before_change,
        #                                  _pdbx_poly_seq_scheme_auth_seq_num_after_change) if i == j])
        # print("Total number of poly_monomers:", len(_pdbx_poly_seq_scheme_auth_seq_num_after_change))
        # print("Amount of changed poly_monomers: ", len(_pdbx_poly_seq_scheme_auth_seq_num_after_change) - same_len)
        # except KeyError:
        # pass
        # print("No_pdbx_poly_seq_scheme.auth_mon_id")

        os.chdir(default_output_path_to_mmCIF)
        io = MMCIFIO()
        io.set_dict(mmcif_dict)
        io.save(mmCIF_ID[:4] + "_out.cif")
        if gzip_mode == 0:
            gzip_for_output_files(mmCIF_ID[:4] + "_out.cif", gzip_mode)
            os.remove(mmCIF_ID[:4] + "_out.cif")
        os.chdir(currentDirectory)
        # print("Processed:", mmCIF_ID[:4])
        # print("File(s) Processed:", stepik)
