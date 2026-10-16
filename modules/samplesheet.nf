/*
 * Sample sheet parser and validator.
 *
 * The sheet is read and checked on the head node, before any job starts, so a
 * bad sheet fails in seconds with every problem listed at once.
 *
 * Columns (header row required):
 *   sample      required  sample name; rows with the same name are lanes and are concatenated
 *   fastq_1     required  R1 FASTQ (gzipped). A relative path is resolved against the
 *                         folder that holds the sample sheet.
 *   fastq_2     required  R2 FASTQ
 *   control     optional  name of an IgG / negative-control sample in the same sheet
 *   is_control  optional  true | false (default false)
 *   pool        optional  library pool ID for read-count scaling (default: one shared pool)
 *
 * Returns a list of maps, one per sample:
 *   [ meta: [id, control, is_control, pool], fastq_1: [Path...], fastq_2: [Path...] ]
 */

import java.nio.file.Files

def parseSamplesheet(sheet_path) {
    def DEFAULT_POOL = 'default'
    def REQUIRED     = ['sample', 'fastq_1', 'fastq_2']
    def KNOWN        = REQUIRED + ['control', 'is_control', 'pool']

    def sheet = file(sheet_path, checkIfExists: true)
    def errors = []

    // ---- Header ----
    def lines = sheet.readLines().findAll { it.trim() }
    if (!lines) {
        error "Sample sheet ${sheet} is empty."
    }
    def header = lines[0].split(',', -1).collect { it.trim() }
    def missing = REQUIRED.findAll { !(it in header) }
    if (missing) {
        error "Sample sheet ${sheet} is missing required column(s): ${missing.join(', ')}.\n" +
              "Expected header: ${KNOWN.join(',')}"
    }
    def unknown = header.findAll { !(it in KNOWN) }
    if (unknown) {
        log.warn "Sample sheet: ignoring unknown column(s): ${unknown.join(', ')}"
    }

    // ---- Rows ----
    def samples   = [:]          // id -> [meta, fastq_1 list, fastq_2 list, first line number]
    def seenFastq = [:]          // FASTQ path -> line number, to catch a lane listed twice
    lines.drop(1).eachWithIndex { line, i ->
        def lineNo = i + 2
        def cells  = line.split(',', -1).collect { it.trim() }
        if (cells.size() != header.size()) {
            errors << "line ${lineNo}: expected ${header.size()} columns, found ${cells.size()}"
            return
        }
        def row = [header, cells].transpose().collectEntries { [(it[0]): it[1]] }

        def id = row.sample
        if (!id) {
            errors << "line ${lineNo}: 'sample' is empty"
            return
        }
        if (!(id ==~ /[A-Za-z0-9._-]+/)) {
            errors << "line ${lineNo}: sample '${id}' may only contain letters, digits, '.', '_' and '-'"
        }

        // is_control
        def isControlText = (row.is_control ?: 'false').toLowerCase()
        if (!(isControlText in ['true', 'false'])) {
            errors << "line ${lineNo}: is_control must be 'true' or 'false', found '${row.is_control}'"
        }
        def meta = [
            id        : id,
            control   : row.control ?: null,
            is_control: isControlText == 'true',
            pool      : row.pool ?: DEFAULT_POOL
        ]

        // FASTQs: each lane needs exactly one R1 and one R2
        def r1 = row.fastq_1
        def r2 = row.fastq_2
        if (!r1 || !r2) {
            errors << "line ${lineNo}: sample '${id}' needs both fastq_1 and fastq_2 " +
                      "(R1/R2 count mismatch: found R1=${r1 ? 1 : 0}, R2=${r2 ? 1 : 0})"
            return
        }
        if (r1 == r2) {
            errors << "line ${lineNo}: fastq_1 and fastq_2 are the same file (${r1})"
        }
        def paths = [r1, r2].collect { resolveFastq(it, sheet) }
        paths.unique(false).each { p ->
            def key = p.toString()
            if (seenFastq.containsKey(key)) {
                errors << "line ${lineNo}: FASTQ ${key} is already listed on line ${seenFastq[key]}"
            } else {
                seenFastq[key] = lineNo
            }
            def problem = checkReadable(p)
            if (problem) {
                errors << "line ${lineNo}: ${problem}"
            }
        }

        // Rows with the same sample name are lanes; their metadata must agree.
        if (samples.containsKey(id)) {
            def first = samples[id]
            ['control', 'is_control', 'pool'].each { k ->
                if (first.meta[k] != meta[k]) {
                    errors << "line ${lineNo}: sample '${id}' has ${k}='${meta[k] ?: ''}', " +
                              "but line ${first.line} has ${k}='${first.meta[k] ?: ''}'"
                }
            }
            first.fastq_1 << paths[0]
            first.fastq_2 << paths[1]
        } else {
            samples[id] = [meta: meta, fastq_1: [paths[0]], fastq_2: [paths[1]], line: lineNo]
        }
    }

    if (!samples && !errors) {
        errors << "the sheet has a header but no sample rows"
    }

    // ---- Controls ----
    samples.each { id, s ->
        def c = s.meta.control
        if (s.meta.is_control && c) {
            errors << "sample '${id}' is a control (is_control=true) but itself lists control '${c}'"
        }
        if (c) {
            if (c == id) {
                errors << "sample '${id}' lists itself as its control"
            } else if (!samples.containsKey(c)) {
                errors << "sample '${id}' lists control '${c}', which is not a sample in this sheet"
            } else if (!samples[c].meta.is_control) {
                errors << "sample '${id}' lists control '${c}', but '${c}' does not have is_control=true"
            }
        }
    }

    if (errors) {
        error "Sample sheet ${sheet} has ${errors.size()} problem(s):\n" +
              errors.collect { "  - ${it}" }.join('\n')
    }

    if (!samples.values().any { it.meta.is_control }) {
        log.info "Sample sheet: no sample has is_control=true. " +
                 "SEACR will run in numeric-threshold mode (top 0.01) for every target."
    }

    return samples.values().collect { [meta: it.meta, fastq_1: it.fastq_1, fastq_2: it.fastq_2] }
}

// Resolve a FASTQ path. Relative paths are taken relative to the sample sheet's folder.
def resolveFastq(String p, sheet) {
    if (p ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ || p.startsWith('/')) {
        return file(p)
    }
    return file(sheet.parent.resolve(p).toString())
}

// Return an error message if a local FASTQ cannot be read, else null.
// Remote paths (s3://, https://, ...) are checked by Nextflow when they are staged.
def checkReadable(p) {
    if (p.scheme != 'file') {
        return null
    }
    if (!Files.exists(p)) {
        return "FASTQ not found: ${p}"
    }
    if (!Files.isRegularFile(p) || !Files.isReadable(p)) {
        return "FASTQ is not a readable file: ${p}"
    }
    return null
}

// Print a short summary table of the parsed sheet.
def samplesheetSummary(samples) {
    def rows = samples.collect { s ->
        String.format('  %-24s %-8s %-20s %-10s %d',
            s.meta.id, s.meta.is_control, s.meta.control ?: '-', s.meta.pool, s.fastq_1.size())
    }
    return String.format('  %-24s %-8s %-20s %-10s %s', 'sample', 'control?', 'control', 'pool', 'lanes') +
           '\n' + rows.join('\n')
}
