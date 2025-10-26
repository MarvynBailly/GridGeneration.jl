# Running the GridGeneration GUI - Quick Reference

## If Nothing Happens When You Launch

The GUI might be running but the window isn't appearing. Here's what to check:

### Step 1: Test GLMakie First
```julia
cd("examples/gui")
include("test_glmakie.jl")
```

**Expected**: A window with a sine wave should appear
**If it works**: GLMakie is fine, proceed to Step 2
**If it doesn't**: See "GLMakie Not Working" below

### Step 2: Run the GUI
```julia
include("GridGenerationGUI.jl")
```

**Watch the terminal for messages**:
```
[ Info: Launching GridGeneration GUI...
[ Info: Initializing GUI components...
[ Info: Building control panel...
[ Info: Opening GUI window...
[ Info: GUI ready! Window should now be visible.
```

**If you see these messages**: The GUI is running!
- Check your taskbar for a new window
- Check if window appeared behind other windows
- Try Alt+Tab to find the window

**If script exits immediately**: The window closed too fast
- This is now fixed with the `wait()` call
- Update your copy if still seeing this

### Step 3: Keep It Running

The GUI needs Julia to keep running to stay open.

**Correct**: Terminal shows "Press Ctrl+C to exit" and waits
**Wrong**: Terminal returns to prompt immediately

To close: Press Ctrl+C in the terminal or close the window

## GLMakie Not Working

### Check 1: Is GLMakie Installed?
```julia
using Pkg
Pkg.status("GLMakie")
```

**Should show**: `GLMakie v0.x.x`
**If not**: Run `Pkg.add("GLMakie")`

### Check 2: OpenGL Support
```julia
using GLMakie
GLMakie.activate!()
```

**No errors**: OpenGL is available
**Errors**: Your system might not support OpenGL

### Check 3: Graphics Drivers
- Update your graphics drivers
- Restart computer after update
- Try again

### Alternative: Use Web GUI Instead
```julia
include("GridGenerationWebGUI.jl")
```

The Web GUI doesn't need OpenGL and might work better on your system.

## Common Issues

### "Window opens then immediately closes"
**Fixed in current version** - we added `wait()` to keep it open

### "Window is frozen/unresponsive"
- The GUI might be loading - wait a few seconds
- Check terminal for errors
- Try resizing the window
- Restart if truly frozen

### "Can't see any controls"
- Window might be too small
- Resize to at least 1400x800
- Default is 1800x1000

### "Nothing happens when I click buttons"
- Domain must be loaded first
- Watch terminal for error messages
- Check status label at bottom of controls

## Platform-Specific

### Windows
- Should work out of the box with recent Windows 10/11
- Update graphics drivers if issues

### macOS
- May need to grant accessibility permissions
- Check System Preferences > Security & Privacy

### Linux
- Needs X11 or Wayland
- Set DISPLAY environment variable if needed:
  ```bash
  export DISPLAY=:0
  ```

## Verification Checklist

Before reporting an issue, verify:

- [ ] GLMakie is installed: `using GLMakie` works
- [ ] test_glmakie.jl shows a window
- [ ] Running from correct directory (examples/gui/)
- [ ] Terminal shows "GUI ready!" message
- [ ] Checked taskbar for window
- [ ] Tried Alt+Tab to find window
- [ ] No errors in terminal output

## Still Not Working?

Try Web GUI instead:
```julia
include("GridGenerationWebGUI.jl")
```

Or use the command-line version:
```julia
cd("../..")
include("examples/generalExample.jl")
```

## Getting Help

Include this info when asking for help:

```julia
# Run this and share output:
using Pkg
Pkg.status("GLMakie")
versioninfo()

# Also share:
# - Operating system
# - Graphics card/driver version
# - Any error messages from terminal
# - Output of test_glmakie.jl
```
