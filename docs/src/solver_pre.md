# Pre-process

```@docs
initialize
```

The pre-process solver initializes the simulation that returns solver set, control volumes, interfaces, and current time.
It could be a new simulation or restart of an interrupted one.
- new run: .txt / .cfg / .toml / etc.
- restart: .jld2
