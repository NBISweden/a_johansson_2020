Help on module gwasim:

NAME gwasim

CLASSES builtins.object SettingsSingleton

    class SettingsSingleton(builtins.object)
     |  Singleton class to store user-defined settings.
     |  
     |  Attributes:
     |      bed_regions_path (str): Path to the BED file containing region data.
     |      bgen_file_path (str): Path to the BGEN file.
     |      threshold (float): Alleles with a frequency below this threshold are considered rare.
     |      my_chr (str): Chromosome of interest.
     |      num_common (int): Number of common variant markers used for simulation.
     |      num_rare (int): Number of rare variant markers used for simulation.
     |      rare_eff_mean (float): Mean effect size for rare variants.
     |      rare_eff_std_dev (float): Standard deviation of effect size for rare variants.
     |      common_eff_mean (float): Mean effect size for common variants.
     |      common_eff_std_dev (float): Standard deviation of effect size for common variants.
     |      err_mean (float): Mean for the error term in phenotype simulation.
     |      err_sd (float): Standard deviation for the error term in phenotype simulation.
     |  
     |  Methods defined here:
     |  
     |  init_settings(self)
     |  
     |  ----------------------------------------------------------------------
     |  Static methods defined here:
     |  
     |  __new__(cls)
     |      Create and return a new object.  See help(type) for accurate signature.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)

FUNCTIONS fix_gt(gt) Fix genotype data.

        Args:
            gt (numpy.ndarray): Genotype data.
        
        Returns:
            numpy.ndarray: Fixed genotype data.

    get_maf(variant)
        Calculate Minor Allele Frequency (MAF) from variant data.
        
        Args:
            variant (object): Variant data.
        
        Returns:
            float: MAF value.

    read_bed_file(S)
        Read and parse region data from a BED file.
        
        Args:
            S (SettingsSingleton): Singleton object containing user-defined settings.
        
        Returns:
            list: A list of dictionaries containing region information.

    select_variants(bgen, region, S)
        Select variants for simulation.
        
        Args:
            bgen (PyBGEN): PyBGEN object for reading BGEN data.
            region (dict): Region dictionary containing region information.
            S (SettingsSingleton): Singleton object containing user-defined settings.
        
        Returns:
            dict: A dictionary containing selected variants and metadata for simulation.

    set_settings_from_command_line_args(S)
        Set settings in the SettingsSingleton (S) based on command-line arguments.

    simulate_phenotype(bgen, variants, S)
        Simulate phenotypes based on selected variants.
        
        Args:
            bgen (PyBGEN): PyBGEN object for reading BGEN data.
            variants (dict): Dictionary containing selected variants for simulation.
            S (SettingsSingleton): Singleton object containing user-defined settings.
        
        Returns:
            dict: A dictionary containing simulated phenotype data and metadata.

    validate_regions(regions, S)
        Validate regions based on simulation criteria.
        
        Args:
            regions (list): List of region dictionaries.
            S (SettingsSingleton): Singleton object containing user-defined settings.
        
        Returns:
            list: A filtered list of validated regions.

FILE
/Users/kiero/Dropbox/WABI/Projects/Johansson_PP/a_johansson_2020/workflow/data/gwasim.py