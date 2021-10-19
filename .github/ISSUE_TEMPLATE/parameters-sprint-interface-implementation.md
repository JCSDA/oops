---
name: "Parameters Sprint: Refactoring a class in ufo or model interface repository"
about: Use this template to create a new issue for the sprint
title: "Parameter implementation for (repo/class):"
labels: 'Sprint'
assignees: ''

---

Update the class identified in the title (hereafter referred to as the **target class**) to process yaml configurations with the ``oops::Parameters`` class.

- define a ``oops::Parameters`` subclass for the target class (with documentation!!)
- add a public ``Parameters_`` typedef
- replace ``eckit::Configuration`` objects in target class constructors with equivalent ``oops::Parameters`` objects
- make necessary changes in dependent code that uses these parameters
- run unit tests making use of the target class and verify that they pass. Make sure you can get one of the tests to fail with a YAML validation error by introducing a typo in the name of one of the YAML options included in the new ``Parameters`` subclass.
