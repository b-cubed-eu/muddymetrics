# Tidyverse R Style Guide (CLI Template)

Use this template for generating R code and updating the b3gbi package to ensure consistency with Tidyverse standards.

## 1. Naming Conventions
- **Functions & Variables:** Use `snake_case`. Names should be concise but meaningful.
  - *Example:* `calculate_summary_stats <- function(df_input) { ... }`
- **Avoid Dots:** Do not use `.` in names (e.g., `data.frame`) to avoid S3 method confusion.

## 2. Syntax & Spacing
- **Assignment:** Use `<-`, never `=`.
- **Spacing:**
  - Always put a space after a comma, but never before.
  - Place spaces around all infix operators (`==`, `+`, `-`, `<-`, etc.).
  - *Example:* `x <- f(a, b)`
- **Line Length:** Limit code to 80 characters per line.
- **Indentation:** Use two spaces for indentation. Never use tabs.

## 3. Functions & Pipes
- **The Pipe:** Use the native pipe `|>` or `%>%`.
  - Always have a space before the pipe and a new line after it.
  - Indent the next line by two spaces.
- **Returns:** Use implicit returns for simple functions (omit `return()`). Only use `return()` for early exits.
  - *Example:* `add_one <- function(x) x + 1`
- **Curly Braces:** `{` should be the last character on the line. `}` should be on its own line.

## 4. Documentation (roxygen2)
- Document all functions using `#'` comments.
- Use `@param` for arguments, `@return` for output description, and `@export` if the function is part of the public API.
- Use `@examples` to provide reproducible usage snippets.

## 5. Error Handling
- Use `cli::cli_abort()` for errors and `cli::cli_warn()` for warnings to provide consistent, well-formatted feedback in the console.
