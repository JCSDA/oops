---
name: "Parameters Sprint: Refactoring an oops application"
about: Use this template to create a new issue for the sprint
title: "Parameter implementation for (application):"
labels: 'Sprint'
assignees: ''

---

Update the oops application identified in the title (hereafter referred to as the **target application**) to access YAML configuration options through a subclass of the ``oops::Parameters`` class rather than through an ``eckit::Configuration`` object.

- define a ``oops::Parameters`` subclass containing the application's configuration options (with documentation!)
- create an instance of this subclass at the top of the application's ``execute()`` member function and call ``validateAndDeserialize()`` to validate the input YAML file and load its contents into the ``oops::Parameters`` object
- refactor the rest of the ``execute()`` function to retrieve configuration options from the ``oops::Parameters`` object rather than the input ``eckit::Configuration`` object
- run unit tests executing the target application and verify that they pass. Make sure you can get one of the tests to fail with a YAML validation error by introducing a typo in the name of a YAML option accessed by the application.
