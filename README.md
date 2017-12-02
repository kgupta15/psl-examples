




<h1 align="center">
  <br>
  psl-examples
  <br>
</h1>

<h4 align="center">Probabilistic Soft Logic based examples models for given data.</h4>

<p align="center">

  <a href="https://codeclimate.com/github/daemonslayer/psl-examples">
    <img src="https://codeclimate.com/github/daemonslayer/psl-examples/badges/gpa.svg" alt="code climate">
  </a>
  <a href="https://codeclimate.com/github/daemonslayer/psl-examples/coverage">
    <img src="https://codeclimate.com/github/daemonslayer/psl-examples/badges/coverage.svg" alt="test coverage">
  </a>
  <a href="https://codeclimate.com/github/daemonslayer/psl-examples">
    <img src="https://codeclimate.com/github/daemonslayer/psl-examples/badges/issue_count.svg" alt="issue count">
  </a>  
  <a>
      <img src="https://img.shields.io/github/license/mashape/apistatus.svg" alt="License">
  </a>
</p>

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
    


### License
MIT License

