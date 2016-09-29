##Install JuliaLang
```shell
sudo apt-get install julia
```

##Install Package in Julia
```julia
Pkg.add("DataFrames")
```

###Install script
```shell
nano install.jl
```

```
#!/usr/bin/julia

Pkg.add("DataFrames")
```

```shell
chown a+r install.jl
./install.jl
```
