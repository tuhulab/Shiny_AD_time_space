# Contributing to Shiny AD Time and Space

Thank you for your interest in contributing to this project! This document provides guidelines for contributing.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue with:
- A clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Screenshots (if applicable)
- Your R version and package versions (`sessionInfo()`)

### Suggesting Enhancements

For feature requests or enhancements:
- Use a clear, descriptive title
- Provide a detailed description of the proposed feature
- Explain why this enhancement would be useful
- Include mockups or examples if applicable

### Code Contributions

#### Setting Up Your Development Environment

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR-USERNAME/Shiny_AD_time_space.git
   cd Shiny_AD_time_space
   ```
3. Install dependencies:
   ```r
   source("install_dependencies.R")
   ```
4. Create a branch for your work:
   ```bash
   git checkout -b feature/your-feature-name
   ```

#### Making Changes

1. **Follow the Code Style Guide**: Read [CODE_STYLE.md](CODE_STYLE.md)
2. **Write Clear Commit Messages**: 
   - Use present tense ("Add feature" not "Added feature")
   - Be descriptive but concise
   - Reference issues when applicable
3. **Test Your Changes**:
   ```r
   source("test_app.R")
   ```
4. **Update Documentation**:
   - Update README.md if needed
   - Add comments to complex code
   - Update CHANGELOG.md

#### Code Quality

Before submitting:

1. Run the code quality check:
   ```r
   source("check_code_quality.R")
   ```
2. Ensure no new linting issues are introduced
3. Test the application locally:
   ```r
   shiny::runApp("App")
   ```

#### Pull Request Process

1. Update CHANGELOG.md with your changes
2. Ensure all tests pass
3. Push to your fork
4. Create a Pull Request with:
   - Clear title and description
   - Reference to any related issues
   - Screenshots of UI changes (if applicable)
   - Summary of changes made

### Code Review

All submissions require review. We use GitHub pull requests for this purpose.

## Development Guidelines

### Project Structure

```
Shiny_AD_time_space/
├── App/
│   ├── ad_in_time_space.R  # Main application
│   ├── config.R            # Configuration
│   └── helpers.R           # Helper functions
├── data/                   # Data files
├── DESCRIPTION             # Package dependencies
├── README.md               # Project documentation
└── test_app.R              # Testing script
```

### Adding New Features

When adding new features:

1. **Configuration**: Add constants to `App/config.R`
2. **Helper Functions**: Add reusable logic to `App/helpers.R`
3. **UI Changes**: Update UI section in `App/ad_in_time_space.R`
4. **Server Logic**: Update server section in `App/ad_in_time_space.R`
5. **Documentation**: Update relevant documentation files

### Adding New Dependencies

If you need to add a new package:

1. Add it to `DESCRIPTION`
2. Add it to `install_dependencies.R`
3. Document why it's needed in your PR

### Testing

While we don't have formal unit tests yet, please:

1. Test your changes manually
2. Run `test_app.R` to verify basic functionality
3. Test with different inputs and edge cases
4. Verify data loading works correctly

## Style Guidelines

### R Code Style

- Use 2 spaces for indentation
- Maximum line length: 80 characters (when practical)
- Use `snake_case` for variables and functions
- Use `SCREAMING_SNAKE_CASE` for constants
- Add spaces around operators
- Follow tidyverse style guide

### Documentation Style

- Use Roxygen-style comments for functions
- Add section headers with `# ----`
- Comment complex logic
- Keep comments up-to-date

### Commit Message Style

```
type: brief description

Longer explanation if needed.
Can span multiple lines.

Fixes #123
```

Types: feat, fix, docs, style, refactor, test, chore

## Questions?

If you have questions, please:
1. Check existing issues
2. Review documentation
3. Contact: Tu Hu (UYHDK AT leo-pharma DOT com)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
