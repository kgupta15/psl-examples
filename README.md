# psl-model
Probabilistic Soft Logic based models for given data

### How to run a PSL based Model

1. Change to the top-level directory of its project (the directory with the Maven `pom.xml` file).

2. Compile the project

    ```
    mvn compile
    ```
    
3. Now use Maven to generate a classpath for your project's dependencies

    ```
    mvn dependency:build-classpath -Dmdep.outputFile=classpath.out
    ```

4. You can now run a class with the command

    ```
    java -cp ./target/classes:`cat classpath.out` <fully qualified class name>
    ```
    
