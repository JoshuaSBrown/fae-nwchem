"""
NWChem Input File Parser

A Python module to parse NWChem input files (.nw) and extract calculation
parameters including basis set, method, and other settings.

This parser handles the standard NWChem input file format with directives
like basis, task, dft, scf, mp2, ccsd, tce, geometry, etc.
"""

import re
from pathlib import Path
from typing import Optional, Union


class NWChemInputParser:
    """
    Parser for NWChem input files (.nw).
    
    Extracts basis set information, calculation method, and other settings
    from NWChem input files.
    
    Example usage:
        parser = NWChemInputParser()
        result = parser.parse("water.nw")
        print(result)
    """
    
    # Supported calculation methods/theories in NWChem
    SUPPORTED_METHODS = [
        'scf', 'dft', 'mp2', 'ccsd', 'ccsd(t)', 'ccsdt', 'ccsdtq',
        'mcscf', 'selci', 'tce', 'tddft', 'sodft', 'pspw', 'band', 'paw',
        'direct_mp2', 'rimp2', 'qmmm', 'md', 'pspw', 'tce'
    ]
    
    # Common DFT functionals
    DFT_FUNCTIONALS = [
        # LDA
        'slater', 'vwn', 'vwn_1', 'vwn_1_rpa', 'vwn_2', 'vwn_3', 'vwn_4', 
        'vwn_5', 'pw91lda',
        # GGA Exchange
        'becke88', 'xpbe96', 'xpw91', 'optx', 'revpbe', 'rpbe',
        # GGA Correlation
        'perdew86', 'lyp', 'pw91', 'cpbe96', 'op', 'hcth', 'hcth120', 
        'hcth147', 'hcth407', 'hcthp14', 'bc95',
        # Hybrid
        'b3lyp', 'pbe0', 'bhlyp', 'bhandh', 'bhandhlyp', 'acm',
        # Meta-GGA
        'm05', 'm05-2x', 'm06', 'm06-l', 'm06-2x', 'm06-hf', 'm08-hx', 
        'm08-so', 'm11', 'm11-l', 'tpss', 'scan', 'rscan',
        # Range-separated
        'cam-b3lyp', 'lc-blyp', 'lc-pbe', 'lc-wpbe', 'lc-wpbeh', 'wb97', 
        'wb97x', 'wb97x-d', 'wb97x-d3',
        # Double hybrid
        'b2plyp', 'b2gp-plyp',
        # Dispersion corrected
        'b3lyp-d3', 'pbe0-d3', 'b97-d', 'b97-d3',
    ]
    
    def __init__(self):
        """Initialize the parser."""
        self._content = ""
        self._lines = []
        
    def parse(self, file_path: Union[str, Path]) -> dict:
        """
        Parse an NWChem input file and extract calculation information.
        
        Args:
            file_path: Path to the .nw input file
            
        Returns:
            Dictionary containing parsed information with keys:
                - basis_set: dict with basis set information
                - method: str or dict with calculation method details
                - tasks: list of task directives found
                - title: str title of the calculation (if present)
                - charge: int molecular charge (if specified)
                - multiplicity: int spin multiplicity (if specified)
                - geometry: dict with geometry information
                - additional_settings: dict with other parsed settings
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")
            
        with open(file_path, 'r') as f:
            self._content = f.read()
            
        # Normalize content: handle comments and continuations
        self._content = self._preprocess_content(self._content)
        self._lines = self._content.split('\n')
        
        result = {
            'basis_set': self._parse_basis(),
            'method': self._parse_method(),
            'tasks': self._parse_tasks(),
            'title': self._parse_title(),
            'charge': self._parse_charge(),
            'multiplicity': self._parse_multiplicity(),
            'geometry': self._parse_geometry(),
            'additional_settings': self._parse_additional_settings(),
        }
        
        # Clean up None values for cleaner output
        result = {k: v for k, v in result.items() if v is not None}
        
        return result
    
    def parse_string(self, content: str) -> dict:
        """
        Parse NWChem input from a string.
        
        Args:
            content: String containing NWChem input
            
        Returns:
            Dictionary containing parsed information
        """
        self._content = self._preprocess_content(content)
        self._lines = self._content.split('\n')
        
        result = {
            'basis_set': self._parse_basis(),
            'method': self._parse_method(),
            'tasks': self._parse_tasks(),
            'title': self._parse_title(),
            'charge': self._parse_charge(),
            'multiplicity': self._parse_multiplicity(),
            'geometry': self._parse_geometry(),
            'additional_settings': self._parse_additional_settings(),
        }
        
        result = {k: v for k, v in result.items() if v is not None}
        
        return result
    
    def _preprocess_content(self, content: str) -> str:
        """
        Preprocess the input content.
        
        - Remove comments (lines starting with # or content after #)
        - Handle line continuations (backslash at end of line)
        - Normalize whitespace
        """
        lines = content.split('\n')
        processed = []
        
        for line in lines:
            # Remove comments (NWChem uses # for comments)
            if '#' in line:
                line = line.split('#')[0]
            line = line.strip()
            if line:
                processed.append(line)
                
        return '\n'.join(processed)
    
    def _find_block(self, block_name: str) -> Optional[str]:
        """
        Find and extract a block of input (e.g., 'basis ... end').
        
        Args:
            block_name: Name of the block to find (e.g., 'basis', 'dft')
            
        Returns:
            String content of the block, or None if not found
        """
        pattern = rf'(?i)^\s*{block_name}\b(.*?)^\s*end\b'
        match = re.search(pattern, self._content, re.MULTILINE | re.DOTALL)
        if match:
            return match.group(0)
        return None
    
    def _parse_basis(self) -> Optional[dict]:
        """
        Parse the basis set block.
        
        Returns:
            Dictionary with basis set information:
                - default: default basis set (if using library)
                - per_element: dict of element-specific basis sets
                - spherical: bool whether spherical harmonics are used
                - cartesian: bool whether cartesian functions are used
        """
        basis_block = self._find_block('basis')
        if not basis_block:
            return None
            
        result = {
            'per_element': {},
            'spherical': False,
            'cartesian': False,
        }
        
        # Check for spherical/cartesian keywords
        if re.search(r'\bspherical\b', basis_block, re.IGNORECASE):
            result['spherical'] = True
        if re.search(r'\bcartesian\b', basis_block, re.IGNORECASE):
            result['cartesian'] = True
            
        # Parse library specifications: "element library basis_name"
        library_pattern = r'(\*|\w+)\s+library\s+([^\n;]+)'
        library_matches = re.findall(library_pattern, basis_block, re.IGNORECASE)
        
        default_basis = None
        for element, basis_name in library_matches:
            basis_name = basis_name.strip()
            # Handle optional "segment" or "file" keywords
            basis_name = re.sub(r'\s+(segment|file)\s+\S+', '', basis_name).strip()
            
            if element == '*':
                default_basis = basis_name
            else:
                result['per_element'][element.capitalize()] = basis_name
                
        if default_basis:
            result['default'] = default_basis
        elif result['per_element']:
            # If all elements have the same basis, set it as default
            bases = list(set(result['per_element'].values()))
            if len(bases) == 1:
                result['default'] = bases[0]
                
        # Handle inline basis set definitions (not library)
        # Look for element followed by basis functions
        inline_pattern = r'^(\w+)\s+([spdfghSPDFGH])\s*$'
        if re.search(inline_pattern, basis_block, re.MULTILINE):
            result['custom_basis'] = True
            
        # Remove empty entries
        if not result['per_element']:
            del result['per_element']
        if not result.get('spherical'):
            del result['spherical']
        if not result.get('cartesian'):
            del result['cartesian']
            
        return result if result else None
    
    def _parse_tasks(self) -> list:
        """
        Parse all task directives.
        
        Returns:
            List of tuples (method, task_type) for each task found
        """
        # Pattern: task method [operation]
        # operation can be: energy, gradient, optimize, hessian, frequencies, etc.
        task_pattern = r'(?i)^\s*task\s+(\S+)(?:\s+(\S+))?'
        matches = re.findall(task_pattern, self._content, re.MULTILINE)
        
        tasks = []
        for match in matches:
            method = match[0].lower() if match[0] else None
            operation = match[1].lower() if len(match) > 1 and match[1] else 'energy'
            if method:
                tasks.append({
                    'theory': method,
                    'operation': operation
                })
                
        return tasks if tasks else None
    
    def _parse_method(self) -> Optional[dict]:
        """
        Parse the calculation method and its settings.
        
        Returns:
            Dictionary with method information including:
                - name: primary method name
                - settings: method-specific settings
        """
        result = {}
        
        # First, get method from tasks
        tasks = self._parse_tasks()
        if tasks:
            # Primary method is from the last (or main) task
            result['primary'] = tasks[-1]['theory']
            
        # Parse DFT block
        dft_info = self._parse_dft_block()
        if dft_info:
            result['dft'] = dft_info
            if 'primary' not in result:
                result['primary'] = 'dft'
                
        # Parse SCF block
        scf_info = self._parse_scf_block()
        if scf_info:
            result['scf'] = scf_info
            
        # Parse MP2 block
        mp2_info = self._parse_mp2_block()
        if mp2_info:
            result['mp2'] = mp2_info
            
        # Parse CCSD block
        ccsd_info = self._parse_ccsd_block()
        if ccsd_info:
            result['ccsd'] = ccsd_info
            
        # Parse TCE block (Tensor Contraction Engine)
        tce_info = self._parse_tce_block()
        if tce_info:
            result['tce'] = tce_info
            
        return result if result else None
    
    def _parse_dft_block(self) -> Optional[dict]:
        """Parse DFT-specific settings."""
        dft_block = self._find_block('dft')
        if not dft_block:
            return None
            
        result = {}
        
        # Parse XC functional
        xc_pattern = r'(?i)\bxc\s+([^\n;]+)'
        xc_match = re.search(xc_pattern, dft_block)
        if xc_match:
            result['xc_functional'] = xc_match.group(1).strip()
        else:
            # Check for direct functional specification
            for func in self.DFT_FUNCTIONALS:
                if re.search(rf'\b{re.escape(func)}\b', dft_block, re.IGNORECASE):
                    result['xc_functional'] = func
                    break
                    
        # Parse grid
        grid_pattern = r'(?i)\bgrid\s+([^\n;]+)'
        grid_match = re.search(grid_pattern, dft_block)
        if grid_match:
            result['grid'] = grid_match.group(1).strip()
            
        # Parse multiplicity (mult)
        mult_pattern = r'(?i)\bmult\s+(\d+)'
        mult_match = re.search(mult_pattern, dft_block)
        if mult_match:
            result['multiplicity'] = int(mult_match.group(1))
            
        # Check for open-shell DFT
        if re.search(r'\bodft\b', dft_block, re.IGNORECASE):
            result['open_shell'] = True
            
        # Parse convergence settings
        #conv_pattern = r'(?i)\bconvergence\s+([^\n;]+)'
        # Parse convergence settings
        #conv_pattern = r""" 
        #    ^convergence\s+energy\s+
        #    (?P<value>
        #        [-+]?
        #        (?:\d+\.\d*|\.\d+|\d+)
        #        (?:[eE][-+]?\d+)?
        #    )
        #    $
        #"""

        conv_pattern = r"""
            ^convergence\s+energy\s+
            (?P<value>
                [-+]?
                (?:\d+\.\d*|\.\d+|\d+)
                (?:[eE][-+]?\d+)?
            )
            \s*$
        """
        
        conv_match = re.search(conv_pattern, dft_block, re.VERBOSE | re.MULTILINE)
        
        if conv_match:
            result['convergence'] = float(conv_match.group('value'))
        #conv_match = re.search(conv_pattern, dft_block)
        print("convergence match")
        print(conv_match)
        #if conv_match:
        #    result['convergence'] = conv_match.group(1).strip()
            
        # Parse dispersion
        disp_pattern = r'(?i)\bdisp\s*(?:vdw\s+)?(\d+)?'
        disp_match = re.search(disp_pattern, dft_block)
        if disp_match:
            result['dispersion'] = True
            if disp_match.group(1):
                result['dispersion_version'] = int(disp_match.group(1))
                
        # Parse max iterations
        maxiter_pattern = r'(?i)\bmaxiter\s+(\d+)'
        maxiter_match = re.search(maxiter_pattern, dft_block)
        if maxiter_match:
            result['max_iterations'] = int(maxiter_match.group(1))
            
        return result if result else None
    
    def _parse_scf_block(self) -> Optional[dict]:
        """Parse SCF-specific settings."""
        scf_block = self._find_block('scf')
        if not scf_block:
            return None
            
        result = {}
        
        # Parse SCF type (rhf, uhf, rohf)
        for scf_type in ['rhf', 'uhf', 'rohf']:
            if re.search(rf'\b{scf_type}\b', scf_block, re.IGNORECASE):
                result['type'] = scf_type.upper()
                break
                
        # Parse max iterations
        maxiter_pattern = r'(?i)\bmaxiter\s+(\d+)'
        maxiter_match = re.search(maxiter_pattern, scf_block)
        if maxiter_match:
            result['max_iterations'] = int(maxiter_match.group(1))
            
        # Parse convergence threshold
        thresh_pattern = r'(?i)\bthresh\s+([^\n;]+)'
        thresh_match = re.search(thresh_pattern, scf_block)
        if thresh_match:
            result['threshold'] = thresh_match.group(1).strip()
            
        # Check for direct SCF
        if re.search(r'\bdirect\b', scf_block, re.IGNORECASE):
            result['direct'] = True
            
        return result if result else None
    
    def _parse_mp2_block(self) -> Optional[dict]:
        """Parse MP2-specific settings."""
        mp2_block = self._find_block('mp2')
        if not mp2_block:
            return None
            
        result = {}
        
        # Check for frozen core
        if re.search(r'\bfreeze\b', mp2_block, re.IGNORECASE):
            result['frozen_core'] = True
            # Try to get specific frozen core settings
            freeze_pattern = r'(?i)\bfreeze\s+(?:atomic|core|virtual)?\s*(\d+)?'
            freeze_match = re.search(freeze_pattern, mp2_block)
            if freeze_match and freeze_match.group(1):
                result['frozen_orbitals'] = int(freeze_match.group(1))
                
        # Check for tight settings
        if re.search(r'\btight\b', mp2_block, re.IGNORECASE):
            result['tight'] = True
            
        # Parse SCS-MP2
        if re.search(r'\bscs\b', mp2_block, re.IGNORECASE):
            result['scs'] = True
            
        return result if result else None
    
    def _parse_ccsd_block(self) -> Optional[dict]:
        """Parse CCSD-specific settings."""
        ccsd_block = self._find_block('ccsd')
        if not ccsd_block:
            return None
            
        result = {}
        
        # Check for frozen core
        if re.search(r'\bfreeze\b', ccsd_block, re.IGNORECASE):
            result['frozen_core'] = True
            
        # Check for perturbative triples
        if re.search(r'\b\(t\)\b', ccsd_block, re.IGNORECASE):
            result['perturbative_triples'] = True
            
        # Parse max iterations
        maxiter_pattern = r'(?i)\bmaxiter\s+(\d+)'
        maxiter_match = re.search(maxiter_pattern, ccsd_block)
        if maxiter_match:
            result['max_iterations'] = int(maxiter_match.group(1))
            
        return result if result else None
    
    def _parse_tce_block(self) -> Optional[dict]:
        """Parse TCE (Tensor Contraction Engine) settings."""
        tce_block = self._find_block('tce')
        if not tce_block:
            return None
            
        result = {}
        
        # Parse method within TCE
        method_pattern = r'(?i)\b(ccsd|ccsd\(t\)|ccsdt|ccsdtq|cisd|cisdt|mbpt2|mbpt3|mbpt4|cr-ccl|eomccsd)\b'
        method_match = re.search(method_pattern, tce_block)
        if method_match:
            result['method'] = method_match.group(1).upper()
            
        # Check for frozen core
        if re.search(r'\bfreeze\b', tce_block, re.IGNORECASE):
            result['frozen_core'] = True
            
        return result if result else None
    
    def _parse_title(self) -> Optional[str]:
        """Parse the title directive."""
        pattern = r'(?i)^\s*title\s+"([^"]+)"'
        match = re.search(pattern, self._content, re.MULTILINE)
        if match:
            return match.group(1)
            
        # Also try without quotes
        pattern = r'(?i)^\s*title\s+([^\n]+)'
        match = re.search(pattern, self._content, re.MULTILINE)
        if match:
            return match.group(1).strip().strip('"\'')
            
        return None
    
    def _parse_charge(self) -> Optional[int]:
        """Parse the molecular charge."""
        # Can be in geometry block or as standalone
        pattern = r'(?i)\bcharge\s+(-?\d+)'
        match = re.search(pattern, self._content)
        if match:
            return int(match.group(1))
        return None
    
    def _parse_multiplicity(self) -> Optional[int]:
        """Parse the spin multiplicity."""
        # Check various places where multiplicity can be specified
        patterns = [
            r'(?i)\bmultiplicity\s+(\d+)',
            r'(?i)\bmult\s+(\d+)',
        ]
        for pattern in patterns:
            match = re.search(pattern, self._content)
            if match:
                return int(match.group(1))
        return None
    
    def _parse_geometry(self) -> Optional[dict]:
        """
        Parse the geometry block.
        
        Returns basic geometry information without full coordinate parsing.
        """
        geom_block = self._find_block('geometry')
        if not geom_block:
            return None
            
        result = {}
        
        # Check for units
        if re.search(r'\bangstroms?\b', geom_block, re.IGNORECASE):
            result['units'] = 'angstrom'
        elif re.search(r'\bbohr\b', geom_block, re.IGNORECASE):
            result['units'] = 'bohr'
        elif re.search(r'\bau\b', geom_block, re.IGNORECASE):
            result['units'] = 'bohr'
        else:
            result['units'] = 'angstrom'  # default
            
        # Check for symmetry
        sym_pattern = r'(?i)\bsymmetry\s+(\S+)'
        sym_match = re.search(sym_pattern, geom_block)
        if sym_match:
            result['symmetry'] = sym_match.group(1)
            
        # Check for noautosym/noautoz
        if re.search(r'\bnoautosym\b', geom_block, re.IGNORECASE):
            result['autosym'] = False
        if re.search(r'\bnoautoz\b', geom_block, re.IGNORECASE):
            result['autoz'] = False
            
        # Count atoms (simple heuristic)
        # Atoms are lines with element symbol followed by coordinates
        atom_pattern = r'^\s*([A-Za-z]{1,2})\s+[-+]?\d*\.?\d+'
        atoms = re.findall(atom_pattern, geom_block, re.MULTILINE)
        if atoms:
            result['num_atoms'] = len(atoms)
            result['elements'] = list(set(atoms))
            
        return result if result else None
    
    def _parse_additional_settings(self) -> Optional[dict]:
        """Parse additional settings and directives."""
        result = {}
        
        # Parse memory
        mem_pattern = r'(?i)^\s*memory\s+(\S+)\s*(\S+)?'
        mem_match = re.search(mem_pattern, self._content, re.MULTILINE)
        if mem_match:
            result['memory'] = {
                'value': mem_match.group(1),
                'unit': mem_match.group(2) if mem_match.group(2) else 'bytes'
            }
            
        # Parse scratch/permanent directories
        scratch_pattern = r'(?i)^\s*scratch_dir\s+(\S+)'
        scratch_match = re.search(scratch_pattern, self._content, re.MULTILINE)
        if scratch_match:
            result['scratch_dir'] = scratch_match.group(1)
            
        perm_pattern = r'(?i)^\s*permanent_dir\s+(\S+)'
        perm_match = re.search(perm_pattern, self._content, re.MULTILINE)
        if perm_match:
            result['permanent_dir'] = perm_match.group(1)
            
        # Parse start/restart
        if re.search(r'(?i)^\s*start\s+', self._content, re.MULTILINE):
            result['job_type'] = 'start'
        elif re.search(r'(?i)^\s*restart\s+', self._content, re.MULTILINE):
            result['job_type'] = 'restart'
            
        # Parse print settings
        print_pattern = r'(?i)^\s*print\s+([^\n]+)'
        print_match = re.search(print_pattern, self._content, re.MULTILINE)
        if print_match:
            result['print'] = print_match.group(1).strip()
            
        # Check for relativistic settings
        if re.search(r'(?i)\brelativistic\b', self._content):
            result['relativistic'] = True
        if re.search(r'(?i)\bzora\b', self._content):
            result['relativistic'] = 'zora'
        if re.search(r'(?i)\bdkh\b', self._content):
            result['relativistic'] = 'dkh'
            
        # Check for solvation (COSMO)
        if self._find_block('cosmo'):
            result['solvation'] = 'cosmo'
            
        return result if result else None


def parse_nwchem_input(file_path: Union[str, Path]) -> dict:
    """
    Convenience function to parse an NWChem input file.
    
    Args:
        file_path: Path to the .nw input file
        
    Returns:
        Dictionary containing parsed information
    """
    parser = NWChemInputParser()
    return parser.parse(file_path)


def parse_nwchem_string(content: str) -> dict:
    """
    Convenience function to parse NWChem input from a string.
    
    Args:
        content: String containing NWChem input
        
    Returns:
        Dictionary containing parsed information
    """
    parser = NWChemInputParser()
    return parser.parse_string(content)


# Example usage and demonstration
if __name__ == "__main__":
    import json
    import sys
    
    # Example NWChem input for testing
    example_input = '''
title "Water molecule DFT calculation"

start water

memory 1 gb

geometry units angstrom
  O  0.00000  0.00000  0.11726
  H  0.00000  0.75698 -0.46906
  H  0.00000 -0.75698 -0.46906
end

basis
  * library 6-311++G**
end

dft
  xc b3lyp
  mult 1
  maxiter 100
  grid fine
  disp vdw 3
end

task dft optimize
task dft freq
'''
    
    # If a file path is provided as argument, parse that file
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        try:
            result = parse_nwchem_input(file_path)
            print(json.dumps(result, indent=2))
        except FileNotFoundError as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # Parse the example input
        print("Example NWChem Input Parser Demo")
        print("=" * 50)
        result = parse_nwchem_string(example_input)
        print(json.dumps(result, indent=2))
        
        print("\n" + "=" * 50)
        print("Quick access to key information:")
        print(f"  Basis set: {result.get('basis_set', {}).get('default', 'Not found')}")
        print(f"  Method: {result.get('method', {}).get('primary', 'Not found')}")
        if 'dft' in result.get('method', {}):
            print(f"  XC Functional: {result['method']['dft'].get('xc_functional', 'Not specified')}")
