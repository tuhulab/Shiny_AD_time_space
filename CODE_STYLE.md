# Code Style Guide

This document outlines the coding conventions used in this project.

## General Principles

1. **Readability**: Code should be easy to read and understand
2. **Consistency**: Follow established patterns throughout the codebase
3. **Modularity**: Separate concerns into different files and functions
4. **Documentation**: Comment complex logic and document functions

## File Organization

### Main Application (`App/ad_in_time_space.R`)
- Header with metadata and purpose
- Library imports
- Source helper files
- UI definition
- Server logic
- App runner

### Configuration (`App/config.R`)
- Centralized constants
- URL endpoints
- Color schemes
- Default parameters

### Helpers (`App/helpers.R`)
- Reusable utility functions
- Data loading functions
- Plot creation functions

## Naming Conventions

### Variables
- Use `snake_case` for variable names
- Use descriptive names: `time_variation_data` not `tvd`
- Boolean variables should be prefixed with `is_` or `has_`

### Functions
- Use `snake_case` for function names
- Use verb-noun pattern: `load_rds_data()`, `create_datatable()`
- Private helper functions can be prefixed with `.`

### Constants
- Use `SCREAMING_SNAKE_CASE` for constants
- Group related constants in lists: `DATA_URLS`, `SKIN_TYPE_COLORS`

## Function Documentation

Use Roxygen-style comments for functions:

```r
#' Brief description of function
#'
#' Longer description if needed
#'
#' @param param_name Description of parameter
#' @return Description of return value
function_name <- function(param_name) {
  # Implementation
}
```

## Code Structure

### Spacing
- Use 2 spaces for indentation (not tabs)
- Add spaces around operators: `x <- 5`, not `x<-5`
- Add space after commas: `function(a, b)`, not `function(a,b)`

### Line Length
- Keep lines under 80 characters when possible
- Break long function calls across multiple lines:

```r
result <- some_function(
  param1 = value1,
  param2 = value2,
  param3 = value3
)
```

### Pipes
- Use `%>%` pipe for data transformation chains
- Put each step on a new line
- Indent continued pipes by 2 spaces

```r
data %>%
  filter(condition) %>%
  mutate(new_col = calculation) %>%
  select(relevant_cols)
```

## Error Handling

Always include error handling for:
- Data loading from URLs
- User inputs
- File operations

Use `tryCatch()` with informative error messages:

```r
tryCatch(
  risky_operation(),
  error = function(e) {
    message(sprintf("Error: %s", e$message))
    NULL
  }
)
```

## Comments

### Inline Comments
- Use `#` for inline comments
- Place above the code they describe
- Keep comments up-to-date with code changes

### Section Headers
- Use decorative comments for major sections:

```r
# ==============================================================================
# Section Name
# ==============================================================================
```

- Use `# ----` for subsections in RStudio (creates outline)

```r
# Data Loading ----
```

## Shiny-Specific Guidelines

### Reactive Programming
- Keep reactive expressions simple and focused
- Use descriptive names for reactive values
- Document reactive dependencies

### UI/Server Separation
- Keep UI definition separate from server logic
- Use helper functions to reduce UI code complexity
- Modularize complex UI sections

### Performance
- Cache data when appropriate
- Use `isolate()` to prevent unnecessary reactivity
- Consider `reactiveVal()` for simple state management

## Testing

- Test helper functions independently
- Verify data loading functions handle errors
- Test UI components with different inputs
- Document expected behavior in comments

## Version Control

### Commits
- Write clear, descriptive commit messages
- Use present tense: "Add feature" not "Added feature"
- Reference issues when applicable

### Branches
- Use descriptive branch names: `feature/add-cell-type-plot`
- Keep branches focused on single features/fixes

## Dependencies

- Document all package dependencies in DESCRIPTION file
- Use specific package versions when reproducibility is critical
- Consider using `renv` for project-level package management

## Documentation

### README
- Keep updated with current functionality
- Include setup instructions
- Provide usage examples
- Document known issues

### Code Documentation
- Document non-obvious logic
- Explain "why" not just "what"
- Include references to relevant research/papers
