---
title: "1. Getting started"
teaching: 15
exercises: 15
questions:
- How can I get started with metabolic models?
objectives:
- Learn the intuition behind flux balance analysis using a toy model
- Learn how to use FBA to design cell factories using the ABC model
keypoints:
- "`load_model` can be used to load a model."
- "`model.solve()` simulates the model."
- "`reaction.flux` contains the computed flux of the reaction."
- "`model.solve()` also returns a solution object that contains the simulation results."
---

{% include nbviewer_iframe.html %}
