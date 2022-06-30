---
name: Leads Approval Checklist
about: Use this template for substantial code changes that may affect all JEDI Applications
title: "[Checklist]"
labels: 'needs review'
assignees: ['rhoneyager','guillaumevernieres','travissluka','amfox37','ss421','svahl991','jjguerrette','byoung-joo','nancylbaker','sking112']

---

## Announcement

We are planning a substantial change in the JEDI code and/or the jedi stack that will likely affect application workloads and/or automated testing (see below for details).

With this checklist, we are asking that JCSDA project leads, JCSDA application co-owners, and AOP task leads confirm that their teams are ready for this change.

When we have confirmed that all teams are ready, we will make the change.

## Deadline

dd / mm / yy

Please either confirm that your team is ready for the change by the deadline above or request an extension.  If we do not hear from you by this date we will assume that your team is ready for the change.


## Description of change

> enter short description and reference associated PRs
> It may be one or more code changes or a change in infrastructure (e.g. modules, CI)

Associated Pull requests:

- [ ] PR 1
- [ ] ...

## Approval

JCSDA project leads, JCSDA Application co-owners, and AOP task leads:

Please check the box next to your application when your team has tested the change on the platform(s) where you run your workflows.  If your application or experiment is not in the list, please add it.

- [ ] JEDI-GDAS
- [ ] JEDI-GEOS
- [ ] JEDI-UFS
- [ ] SOCA / JEDI-GODAS
- [ ] LAND
- [ ] UM-JEDI
- [ ] LFRic-JEDI
- [ ] MPAS-JEDI
- [ ] Neptune-JEDI
- [ ] Shallow-water

## Systems tested

Format for entries (add entries as needed) is `bundle, modules, date, name` where `name` is the name of the person who tested it (we encourage multiple entries with different testers for each platform/bundle/module combination).  Provide details in the comments below if needed.

Hera:

- [ ] fv3-bundle, jedi/intel-impi, dd/mm, `<name>`
- [ ] soca-science, Hera, jedi/intel-mpi, dd/mm, `<name>`
- [ ] ...

Orion:

- [ ] fv3-bundle, jedi/intel-impi, dd/mm, `<name>`
- [ ] soca-science, Hera, jedi/intel-mpi, dd/mm, `<name>`
- [ ] ...

Discover:

- [ ] fv3-bundle, jedi/intel-impi, dd/mm, `<name>`
- [ ] soca-science, Hera, jedi/intel-mpi, dd/mm, `<name>`
- [ ] ...

S4:

- [ ] fv3-bundle, jedi/intel-impi, dd/mm, `<name>`
- [ ] soca-science, Hera, jedi/intel-mpi, dd/mm, `<name>`
- [ ] ...

Mac:

- [ ] fv3-bundle, clang/openmpi, dd/mm, `<name>`
- [ ] soca-science, Hera, clang/openmpi, dd/mm, `<name>`
- [ ] ...

Other:

- [ ] `<platform>` fv3-bundle, jedi/intel-impi, dd/mm, `<name>`
- [ ] `<platform>` soca-science, jedi/intel-mpi, dd/mm, `<name>`
- [ ] ...

## Rollback plan

> how can we undo this change if it adversely affects our workflow?