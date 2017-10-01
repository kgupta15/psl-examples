## Setup

First, ensure that your system meets the [[prerequisites | Prerequisites]].
Then clone the `psl-examples` repository:
```
git clone https://github.com/daemonslayer/psl-examples.git
```

## Running

Then move into the root directory for the easy link prediction example:
```
cd psl-examples/link_prediction/groovy
```

Each example comes with a `run.sh` script to quickly compile and run the example.
To compile and run the example:
```
./run.sh
```

To see the output of the example, check the `output/default/knows_infer.txt` file:
```
cat output/default/knows_infer.txt
```

You should see some output like:
```
--- Atoms:
KNOWS(Steve, Ben) Truth=[0.64]
KNOWS(Alex, Dhanya) Truth=[0.36]
< ... 48 rows omitted for brevity ...>
KNOWS(Dhanya, Ben) Truth=[0.55]
KNOWS(Alex, Sabina) Truth=[0.44]
# Atoms: 52
```
The exact order of the output may change and some rows were left out for brevity.

Now that we have the example running, lets take a look inside the only source file for the example:
`src/main/java/edu/ucsc/linqs/psl/example/easylp/EasyLP.groovy`.

## Configuration

One of the first things you may notice in this file is a private classes that holds configuration data.
In addition to ConfigBundles, it is sometimes also useful to create configuration classes that you can pass quickly change and run different experiments with.
This is not required, but you may find it useful.

In the `populateConfigBundle()` method you can see the ConfigBundle actually getting created:
```groovy
ConfigBundle cb = ConfigManager.getManager().getBundle("easylp");
```

## Defining Predicates
The `definePredicates()` method defines the three predicates for our example:
```groovy
model.add predicate: "Lived", types: [ConstantType.UniqueID, ConstantType.UniqueID];
model.add predicate: "Likes", types: [ConstantType.UniqueID, ConstantType.UniqueID];
model.add predicate: "Knows", types: [ConstantType.UniqueID, ConstantType.UniqueID];
```

Each predicate here takes two unique identifiers as arguments.
- **Lived** indicates that a person has lived in a specific location. For example: Lived(Sammy, SantaCruz) would indicate that Sammy has lived in Santa Cruz.
- **Likes** indicates the extent to which a person likes something. For example: Likes(Sammy, Hiking) would indicate the extent that Sammy likes hiking.
- **Knows** indicates that a person knows some other person. For example: Knows(Sammy, Jay) would indicate that Sammy and Jay know each other.

## Defining Rules
The `defineRules()` method defines seven rules for the example.
There are pages that cover the PSL [[rule specification | Rule Specification]] and the [[rule specification in Groovy | Rule Specification in Groovy]].
We will discuss the following two rules:
```groovy
model.add(
   rule: ( Lived(P1,L) & Lived(P2,L) & (P1-P2) ) >> Knows(P1,P2),
   squared: config.sqPotentials,
   weight : config.weightMap["Lived"]
);

model.add(
   rule: ~Knows(P1,P2),
   squared:config.sqPotentials,
   weight: config.weightMap["Prior"]
);
```

The first first rule can be read as "If P1 and P2 are different people and have both lived in the same location, L, then they know each other".
Some key points to note from this rule are:
 - The variable `L` was reused in both `Lived` atoms and therefore must refer to the same location.
 - `(P1 - P2)` is shorthand for P1 and P2 referring to different people (different unique ids).

The second rule is a special rule that acts as a prior.
Notice how this rule is not an implication like all the other rules.
Instead, this rule can be read as "By default, people do not know each other".
Therefore, the program will start with the belief that no one knows each other and this prior belief will be overcome with evidence.

## Loading Data
The `loadData()` method loads the data from the flat files in the `data` directory into the data store that PSL is working with.
For berevity, we will only be looking at two files:
```groovy
Inserter inserter = ds.getInserter(Lived, obsPartition);
InserterUtils.loadDelimitedData(inserter, Paths.get(config.dataPath, "lived_obs.txt").toString());

inserter = ds.getInserter(Likes, obsPartition);
InserterUtils.loadDelimitedDataTruth(inserter, Paths.get(config.dataPath, "likes_obs.txt").toString());
```

Both portions load data using the `InserterUtils`.
The primary difference between the two calls is that the second one is looking for a truth value while the first one assumes that 1 is the truth value.

If we look in the files, we see lines like:

`data/lives_obs.txt`
```
Jay	Maryland
Jay	California
```

`data/likes_obs.txt`
```
Jay	Machine Learning  1
Jay	Skeeball 0.8
```

In `lives_obs.txt`, there is no need to use a truth value because living somewhere is a discrete act.
You have either lived there or you have not.
Liking something, however, is more continuous.
Jay may like *Machine Learning* 100%, but he only likes *Skeeball* 80%.

### Partitions
Here we must take a moment to talk about data partitions.
In PSL, we use **partitions** to organize data.
A partition is nothing more than a container for data, but we use them to keep specific chunks of data together or separate.
For example if we are running evaluation, we must be sure not use our test partition in training.

PSL users typically organize their data in at least three different partitions (all of which you can see in this example):
- **observations** (called `obsPartition` in this example): In this partition we put actual observed data. In this example, we put all the observations about who has lived where, who likes what, and who knows who in the observations partition.
- **targets** (called `targetsPartition` in this example): In this partition we put atoms that we want to infer values for. For example if we want to if Jay and Sammy know each other, then we would put the atom `Knows(Jay, Sammy)` into the targets partition.
- **truth** (called `truthPartition` in this example): In this partition we put our test set, data that we have actual values for but are not including in our observations for the purpose of evaluation. For example, if we know that Jay and Sammy actually do know each other, we would put `Knows(Jay, Sammy)` in the truth partition with a truth value of 1.

## Running Inference
The `runInference()` method handles running inference for all the data we have loaded.

Before we run inference, we have to set up a database to use for inference:
```groovy
HashSet closed = new HashSet<StandardPredicate>([Lived, Likes]);
Database inferDB = ds.getDatabase(targetsPartition, closed, obsPartition);
```
The `getDatabase()` method of `DataStore` is the proper way to get a database.
This method takes a minimum of two parameters:
- The partition that this database is allowed to write to. In inference, we will be writing the inferred truth values of atom to the target partition, so we will need to have it open for writing.
- A set of partitions to be closed in the write partition. Even though we are writing values into the write partition, we may only have a few predicates that we actually want to infer values for. This parameter allows you to close those predicates that you do not what changed.
Lastly, `getDatabase()` takes any number of read-only partitions that you want to include in this database.
In our example, we want to include our observations when we run inference.

Now we are ready to run inference:
```groovy
MPEInference mpe = new MPEInference(model, inferDB, config.cb);
mpe.mpeInference();
mpe.close();
inferDB.close();
```

To the `MPEInference` constructor, we supply our model, the database to infer over, and our ConfigBundle.
To see the results, then we will need to look inside of the target partition.

## Output
The method `writeOutput()` handles printing out the results of the inference.
There are two key lines in this method:
```groovy
Database resultsDB = ds.getDatabase(targetsPartition);
...
Set atomSet = Queries.getAllAtoms(resultsDB, Knows);
```

The first line gets a fresh database that we can get the atoms from.
Notice that we are passing in `targetsPartition` as a write partition, but we are actually just reading from it.

The second line uses the `Queries` class to get all the `Knows` atoms from the database we just created.

## Evaluation
Lastly, the `evalResults()` method handles seeing how well our model did.
The `DiscretePredictionComparator` and `ContinuousPredictionComparator` classes provide basic tools to compare two partitions.
In this example, we are comparing our target partition to our truth partition.
