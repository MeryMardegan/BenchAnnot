// ============================================================
// Generic helpers
// ============================================================

// Generic validation for directory-type parameters.
// Used by any workflow that needs to check --param paths before running.
def validateRequiredDirs(List<String> paramNames) {
    def missing = paramNames.findAll { name -> !params[name] }
    if (missing) {
        error "Missing required parameters: ${missing.collect { name -> "--${name}" }.join(', ')}"
    }
    paramNames.each { name ->
        def location = file(params[name])
        if (!location.exists() || !location.isDirectory()) {
            error "Parameter --${name} must be an existing directory: ${params[name]}"
        }
    }
}

def databaseImageIsValid(String imagePath) {
    if (!imagePath) {
        return false
    }
    def image = file(imagePath)

    return image.exists() &&
           image.isFile() &&
           image.size() > 0 &&
           image.name.endsWith('.sqsh')
}

// ============================================================
// Resolve the path to a samplesheet file.
// If the path is absolute, use it as is. Otherwise, resolve it relative to the project directory.
// ============================================================
def resolveSamplesheetPath(String rawPath) {
    def file = file(rawPath)
    return file.isAbsolute() ? file : file("${projectDir}/${rawPath}")
}

// ============================================================
// Samplesheet row validation
// ============================================================

// Validate one row from the prokaryotic samplesheet.
// Required metadata must be present, genome FASTA must exist,
// taxid must be numeric, and genetic_code defaults to 11 when absent.
def validateProkaryoteRow(Map row) {
    def requiredFields = ['sample_id', 'genome_fasta', 'species', 'taxid']
    def missing = requiredFields.findAll { field -> !row[field] }
    if (missing) {
        error "Prokaryote samplesheet: missing required field(s) [${missing.join(', ')}]. Row: '${row}'"
    }

    def genomeFile = resolveSamplesheetPath(row.genome_fasta)
    if (!genomeFile.exists() || !genomeFile.isFile()) {
        error "Prokaryote samplesheet: genome_fasta not found for sample '${row.sample_id}': '${row.genome_fasta}'"
    }

    def taxid = row.taxid.toString()
    if (!taxid.isInteger()) {
        error "Prokaryote samplesheet: taxid must be numeric for sample '${row.sample_id}': '${row.taxid}'"
    }

    def geneticCode = row.genetic_code
        ? row.genetic_code.toString()
        : '11'
    if (!geneticCode.isInteger()) {
        error "Prokaryote samplesheet: genetic_code must be numeric for sample '${row.sample_id}': '${geneticCode}'"
    }

    return tuple(row.sample_id, genomeFile, row.species, taxid, geneticCode
    )
}

// Validate one row from the eukaryotic samplesheet.
// Required metadata must be present and all referenced files
// must exist on disk.
def validateEukaryoteRow(Map row) {

    def requiredFields = ['sample_id', 'genome_fasta', 'reference_gff', 'reference_faa', 'organism_id']
    def missing = requiredFields.findAll { field -> !row[field] }

    if (missing) {
        error "Eukaryote samplesheet: missing required field(s) [${missing.join(', ')}]. Row:'${row}'"
    }

    def genomeFile = resolveSamplesheetPath(row.genome_fasta)
    def gffFile    = resolveSamplesheetPath(row.reference_gff)
    def faaFile    = resolveSamplesheetPath(row.reference_faa)

    def inputFiles = [
        genome_fasta : genomeFile,
        reference_gff: gffFile,
        reference_faa: faaFile
    ]

    inputFiles.each { label, inputFile ->

        if (!inputFile.exists() || !inputFile.isFile()) {
            error "Eukaryote samplesheet: ${label} not found for sample '${row.sample_id}': ${inputFile}"
        }
    }

    return tuple(
        row.sample_id, genomeFile, gffFile, faaFile, row.organism_id
    )
}

// ============================================================
// Per-tool external resource validation
// ============================================================

def validateBaktaInputs() {
    validateRequiredFiles(['bakta_db'])
}

def validateKofamInputs() {
    validateRequiredDirs(['kofamscan_db'])
}

def validateInterproscanInputs() {
    validateRequiredFiles(['ips_data_dir'])
}

def validatePgapInputs() {
    validateRequiredDirs(['pgap_dir'])
    validateRequiredFiles(['pgap_container'])
}

// ============================================================
// Generic validation for file-type parameters
// ============================================================

def validateRequiredFiles(List<String> paramNames) {

    def missing = paramNames.findAll { name -> !params[name] }

    if (missing) {
        error "Missing required parameters: ${missing.collect { name -> "--${name}" }.join(', ')}"
    }

    paramNames.each { name ->

        def location = file(params[name])

        if (!location.exists() || !location.isFile() || location.size() <= 0) {
            error "Parameter --${name} must be an existing non-empty file: ${params[name]}"
        }
    }
}
