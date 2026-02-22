# Create diff_full.txt
old <- "original_version/rethinking_type_s_and_m_errors.qmd"
new <- "rethinking_type_s_and_m_errors.qmd"
diff_out <- "diff.txt"

word_rx <- "@[-[:alnum:]_:.]+|[[:alnum:]_]+|[^[:space:]]"

args <- c(
  "-c", "core.autocrlf=false",
  "-c", "diff.algorithm=patience",
  "--no-pager", "diff", "--no-index",
  "-w", "--ignore-cr-at-eol",
  "--word-diff=plain",
  paste0("--word-diff-regex=", word_rx),
  "-U100000",
  old, new
)

status <- system2("git", args, stdout = diff_out)
if (!is.null(attr(status, "status")) && status != 0) {
  stop("git diff failed with status ", status)
}
cat("Wrote diff to:", normalizePath(diff_out), "\n")


render_qmd_worddiff_html <- function(input_path,
                                     output_path,
                                     title = "Diff (QMD)",
                                     show_line_numbers = TRUE) {
  
  # ---- read ----------------------------------------------------------------
  raw <- readLines(input_path, warn = FALSE)
  raw <- sub("\r$", "", raw)  # normalize CRLF
  n_raw <- length(raw)
  
  # ---- helpers --------------------------------------------------------------
  strip_tokens <- function(x) gsub("\\{\\+|\\+\\}|\\[\\-|\\-\\]", "", x, perl = TRUE)
  
  html_escape <- function(x) {
    x <- gsub("&", "&amp;",  x, fixed = TRUE)
    x <- gsub("<", "&lt;",   x, fixed = TRUE)
    x <- gsub(">", "&gt;",   x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
  }
  
  # Apply diff highlighting on a WHOLE TEXT blob (handles multi-line markers)
  # FIX: avoid treating Quarto/Pandoc citations like [-@key] as deletions
  apply_worddiff_spans_blob <- function(x) {
    x <- gsub("(?s)\\{\\+(.*?)\\+\\}", "<span class=\"ins\">\\1</span>", x, perl = TRUE)
    x <- gsub("(?s)\\[\\-(?!@)(.*?)\\-\\]", "<span class=\"del\">\\1</span>", x, perl = TRUE)
    x
  }
  
  # Neutralize punctuation-only citation wrapper swaps:
  #   [-(-]{+[+}  -> [
  #   [-)-]{+]+}  -> ]
  # This removes red/green for ( ) vs [ ] around citations like:
  #   ... tests [-([see also @key)].  ->  ... tests [see also @key].
  collapse_cite_bracket_swaps <- function(x) {
    x <- gsub("(?s)<span class=\"del\">\\(</span>\\s*<span class=\"ins\">\\[</span>", "[", x, perl = TRUE)
    x <- gsub("(?s)<span class=\"del\">\\)</span>\\s*<span class=\"ins\">\\]</span>", "]", x, perl = TRUE)
    x
  }
  
  # Conservative inline code outside fenced blocks
  apply_inline_code <- function(x) {
    gsub("`([^`]+)`", "<code class=\"icode\">\\1</code>", x, perl = TRUE)
  }
  
  is_fence_line <- function(clean_line) {
    grepl("^\\s*`{3,}(\\{[^}]*\\})?\\s*$", clean_line, perl = TRUE)
  }
  
  is_diff_meta <- function(clean_line) {
    grepl("^(diff --git|index\\s+|@@\\s|---\\s|\\+\\+\\+\\s)", clean_line, perl = TRUE)
  }
  
  # OPTION B: only headers if the RAW line starts with '#'
  is_real_md_header_raw <- function(raw_line) {
    grepl("^\\s*#{1,6}\\s+\\S", raw_line, perl = TRUE)
  }
  
  parse_header_level_raw <- function(raw_line) {
    m <- regexec("^\\s*(#{1,6})\\s+.*$", raw_line, perl = TRUE)
    hit <- regmatches(raw_line, m)[[1]]
    if (length(hit) >= 2) nchar(hit[2]) else NA_integer_
  }
  
  # ---- pre-process text safely ---------------------------------------------
  blob <- paste(raw, collapse = "\n")
  blob <- html_escape(blob)
  blob <- apply_worddiff_spans_blob(blob)
  blob <- collapse_cite_bracket_swaps(blob)
  
  lines_html <- strsplit(blob, "\n", fixed = TRUE)[[1]]
  
  # Fallback if something odd happened (shouldn't, but humans love edge cases)
  if (length(lines_html) != n_raw) {
    lines_html <- vapply(raw, function(s) {
      s2 <- html_escape(s)
      s2 <- apply_worddiff_spans_blob(s2)
      collapse_cite_bracket_swaps(s2)
    }, character(1))
  }
  
  # ---- render blocks --------------------------------------------------------
  out <- character(0)
  in_code <- FALSE
  code_buf <- character(0)
  
  flush_code <- function() {
    out <<- c(out,
              "<pre class=\"codeblock\">",
              paste(code_buf, collapse = "\n"),
              "</pre>")
    code_buf <<- character(0)
  }
  
  emit_line_div <- function(html_line, klass = "line", n = NULL) {
    if (show_line_numbers && !is.null(n)) {
      paste0(
        "<div class=\"", klass, "\">",
        "<span class=\"ln\">", n, "</span>",
        "<span class=\"txt\">", html_line, "</span>",
        "</div>"
      )
    } else {
      paste0("<div class=\"", klass, "\">", html_line, "</div>")
    }
  }
  
  for (i in seq_len(n_raw)) {
    clean <- strip_tokens(raw[i])
    
    # Fence toggling: use clean line so {+```+} / [-```-] still works
    if (is_fence_line(clean)) {
      if (!in_code) {
        in_code <- TRUE
        code_buf <- character(0)
        code_buf <- c(code_buf, lines_html[i])
      } else {
        code_buf <- c(code_buf, lines_html[i])
        flush_code()
        in_code <- FALSE
      }
      next
    }
    
    if (in_code) {
      code_buf <- c(code_buf, lines_html[i])
      next
    }
    
    # Outside code: diff meta lines
    if (is_diff_meta(clean)) {
      out <- c(out, emit_line_div(lines_html[i], klass = "meta", n = i))
      next
    }
    
    # Outside code: headers (Option B)
    if (is_real_md_header_raw(raw[i])) {
      lvl <- parse_header_level_raw(raw[i])
      
      header_text_raw <- sub("^\\s*#{1,6}\\s+", "", raw[i], perl = TRUE)
      header_text_raw <- sub("^\\*\\*(.*?)\\*\\*$", "\\1", header_text_raw, perl = TRUE)
      header_text_raw <- sub("^_(.*?)_$", "\\1", header_text_raw, perl = TRUE)
      
      header_text_html <- collapse_cite_bracket_swaps(
        apply_worddiff_spans_blob(html_escape(header_text_raw))
      )
      
      out <- c(out, sprintf("<h%d>%s</h%d>", lvl, header_text_html, lvl))
      next
    }
    
    # Normal text line
    line <- apply_inline_code(lines_html[i])
    out <- c(out, emit_line_div(line, klass = "line", n = i))
  }
  
  # Close dangling code block if file ends mid-block
  if (in_code) flush_code()
  
  # ---- final HTML -----------------------------------------------------------
  css <- "
:root{--ins:#eaffea;--del:#ffecec;--codebg:#f7f7f7;--metabg:#f3f3f3;--border:#e8e8e8;}
body{font-family:'Times New Roman', Times, serif;margin:1.2rem;max-width:1100px;}
h1,h2,h3,h4,h5,h6{margin:0.9em 0 0.4em;line-height:1.2;}
h1{font-size:2.0em;} h2{font-size:1.6em;} h3{font-size:1.35em;}
h4{font-size:1.2em;} h5{font-size:1.05em;} h6{font-size:1.0em;}

.line,.meta{white-space:pre-wrap;line-height:1.5;}
.meta{background:var(--metabg);border:1px solid var(--border);padding:6px 8px;border-radius:4px;margin:6px 0;}
.ins{background:var(--ins);}
.del{background:var(--del);text-decoration:line-through;}

pre.codeblock{
  white-space:pre-wrap;line-height:1.35;
  font-family:'Courier New', Courier, monospace;
  background:var(--codebg);border:1px solid var(--border);
  padding:10px 12px;margin:10px 0;border-radius:6px;overflow:auto;
}
code.icode{
  font-family:'Courier New', Courier, monospace;
  background:#f1f1f1;padding:0 3px;border-radius:3px;
}

.ln{
  display:inline-block;min-width:4.5em;
  color:#777;user-select:none;padding-right:0.8em;
}
.txt{display:inline;}
"
html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\">",
  "<title>", html_escape(title), "</title>",
  "<style>", css, "</style></head><body>",
  "<h1>", html_escape(title), "</h1>",
  paste(out, collapse = "\n"),
  "</body></html>"
)

con <- file(output_path, open = "wb")
on.exit(close(con), add = TRUE)
writeBin(charToRaw(enc2utf8(html)), con)

invisible(output_path)
}


render_qmd_worddiff_html(
  input_path  = "diff.txt",
  output_path = "diff_v1_v2.html",
  title       = "QMD word-diff"
)
