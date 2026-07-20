changes made to groot

1. It can read .txt3 files
2. It is capable of doing calibrations


Plugin support

Groot discovers optional plugins without loading their libraries at startup.
Plugin manifests use ROOT TEnv syntax and the `.plugin` extension. Search
directories are read in this order:

1. `GROOT_PLUGIN_PATH` (colon-separated directories)
2. `$HOME/.local/lib/groot/plugins`
3. `$GSYS/lib/groot/plugins`

Example:

  export GROOT_PLUGIN_PATH="$HOME/.local/lib/groot/plugins:/opt/groot/plugins"

A minimal manifest is:

  Plugin.Id: example
  Plugin.Name: Example Plugin
  Plugin.Version: 1.0.0
  Plugin.ApiVersion: 3
  Plugin.Library: libExamplePlugin.so
  Plugin.Action.Count: 1
  Plugin.Action.0.Id: example.run
  Plugin.Action.0.Label: Run example
  Plugin.Action.0.Tooltip: Run the example plugin action

The GUI displays actions from valid manifests. The corresponding library is
loaded only when an action is first selected. Missing libraries, entry symbols,
duplicate identifiers, and incompatible API versions are reported without
terminating Groot.

Plugin API v3 supports one application-overlay session per pad and target.
ROOT native canvas interaction always runs. Active sessions observe normalized
events after ROOT handling while Groot-specific hotkeys, markers, and cursor
policy are suspended for that pad. Sessions close with their owner canvas.
Groot explicitly shuts down plugin sessions and instances before ROOT deletes
remaining canvases, so canvas cleanup remains safe regardless of static object
destruction order.

External Groot connectors compile against the public headers under
`include/Plugin` and the `GPLUGIN` library in Groot's normal build tree.

Running `make` builds Groot without regression-test executables. Run
`make test` to explicitly build and execute the Plugin Manager regression
tests.
