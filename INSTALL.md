# Installation Instructions

## Requirements
  To install this program you need to have **Fortran compiler (gfortran)** and **Make** installed on your system


To run this project locally, follow these steps on a bash-based system:

## 1. Clone the repository:
```bash
git clone https://github.com/stefeversdijk/GoEST
```


## 2. Navigate into the source directory:
```bash
cd EST/src
```

## 3. Compile the program

Use the `Makefile` to compile the program. The `Makefile` provides options for different builds:

  ```bash
  make
  ```

## 5. Return to EST directory
```bash
cd ..
```
## 6. Run the program

```bash
./GoEST "SYSTEM" "BASIS"
```
## Example
./GoEST H2O cc-pvdz
