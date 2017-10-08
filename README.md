# ShapeSearch
Search Similar 3D shapes

NOTE: Right now this program is in very nascent state. I am working on it to make it production ready.

Objectives:
1. This program would be used to compare two or more 3D shapes.
2. At present this works on Linux (tested on ubuntu).
3. I will upload the code for development once I bring it out of embarrassing state.
4. This program would act as tool for people working in CAD/3D printing/Design etc industry to:
    a. Compared shapes.
    b. Search for similar shapes in repository.
    c. Get alerts if you are recreating shapes which already exist in your shape repository.
    
Current State of affairs:
1. Right now it works only on .stl files.
2. No mechnanised error handling.
3. Performance of the tool is extremely bad.
4. Initialization of unsupervised clustering is in bad state.
5. Gives wrong results sometimes.
6. Number of clusters in the shape repository is not calculated automatically.
7. Serialization of formed clusters is missing.
8. Good thing is that if you know the number clusters in your shape repository, it just works!


How to use it:

Using it is quite straight forward:
"./smartModels <<path of shape repository>> <<.stl file to compare with>> <<number of clusters in your repository>>"



