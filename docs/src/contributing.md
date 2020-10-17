# Contributing

This page details the some of the guidelines that should be followed when contributing to this package.

## Branches

Please open pull requests against the `master` branch whenever possible.

### `release-*` branches

  * Each new minor version `x.y.0` gets a branch called `release-x.y` (a [protected branch](https://help.github.com/en/github/administering-a-repository/about-protected-branches)).
  * New versions are usually tagged only from the `release-x.y` branches.

## Style Guide

Follow the style of the surrounding text when making changes. When adding new features please try to stick to the following points whenever applicable.

### Julia

  * 4-space indentation;
  * modules spanning entire files should not be indented, but modules that have surrounding code should;
  * no blank lines at the start or end of files;
  * do not manually align syntax such as `=` or `::` over adjacent lines;
  * use `function ... end` when a method definition contains more than one toplevel expression;
  * related short-form method definitions don't need a new line between them;
  * unrelated or long-form method definitions must have a blank line separating each one;
  * surround all binary operators with whitespace except for `::`, `^`, and `:`;
  * files containing a single `module ... end` must be named after the module;
  * method arguments should be ordered based on the amount of usage within the method body;
  * methods extended from other modules must follow their inherited argument order, not the above rule;
  * explicit `return` should be preferred except in short-form method definitions;
  * avoid dense expressions where possible e.g. prefer nested `if`s over complex nested `?`s;
  * include a trailing `,` in vectors, tuples, or method calls that span several lines;
  * do not use multiline comments (`#=` and `=#`);
  * wrap long lines as near to 92 characters as possible, this includes docstrings;
  * follow the standard naming conventions used in `Base`.

### Markdown

  * Use unbalanced `#` headers, i.e. no `#` on the right hand side of the header text;
  * include a single blank line between toplevel blocks;
  * unordered lists must use `*` bullets with two preceding spaces;
  * do *not* hard wrap lines;
  * use emphasis (`*`) and bold (`**`) sparingly;
  * always use fenced code blocks instead of indented blocks