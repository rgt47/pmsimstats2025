#!/bin/bash
# Navigation Functions Generator
# Creates shell functions for quick directory navigation from anywhere
# Usage: ./navigation_scripts.sh [--install | --uninstall]
#   --install   : Add navigation functions to your shell config
#   --uninstall : Remove navigation functions from your shell config

SHELL_RC="${HOME}/.zshrc"
if [[ "$SHELL" == *"bash"* ]]; then
    SHELL_RC="${HOME}/.bashrc"
fi

# Navigation functions to be added
NAVIGATION_FUNCTIONS='
# ZZCOLLAB Navigation Functions (added by navigation_scripts.sh)
# These allow one-letter navigation from anywhere in your project

# Find project root (looks for DESCRIPTION file)
_zzcollab_root() {
    local dir="$PWD"
    while [[ "$dir" != "/" ]]; do
        if [[ -f "$dir/DESCRIPTION" ]] || [[ -f "$dir/.zzcollab_project" ]]; then
            echo "$dir"
            return 0
        fi
        dir="$(dirname "$dir")"
    done
    return 1
}

# Navigation functions
a() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis" || echo "Not in ZZCOLLAB project"; }
d() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/data" || echo "Not in ZZCOLLAB project"; }
n() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis" || echo "Not in ZZCOLLAB project"; }
f() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis/figures" || echo "Not in ZZCOLLAB project"; }
t() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis/tables" || echo "Not in ZZCOLLAB project"; }
s() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis/scripts" || echo "Not in ZZCOLLAB project"; }
p() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/analysis/paper" || echo "Not in ZZCOLLAB project"; }
r() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root" || echo "Not in ZZCOLLAB project"; }
m() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/man" || echo "Not in ZZCOLLAB project"; }
e() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/tests" || echo "Not in ZZCOLLAB project"; }
o() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/docs" || echo "Not in ZZCOLLAB project"; }
c() { local root=$(_zzcollab_root); [[ -n "$root" ]] && cd "$root/archive" || echo "Not in ZZCOLLAB project"; }

# List navigation shortcuts
nav() {
    echo "ZZCOLLAB Navigation Shortcuts:"
    echo "  r → project root"
    echo "  d → data/"
    echo "  a/n → analysis/"
    echo "  s → analysis/scripts/"
    echo "  p → analysis/paper/"
    echo "  f → analysis/figures/"
    echo "  t → analysis/tables/"
    echo "  m → man/"
    echo "  e → tests/"
    echo "  o → docs/"
    echo "  c → archive/"
}
# End ZZCOLLAB Navigation Functions
'

# Function to install navigation functions
install_functions() {
    if grep -q "ZZCOLLAB Navigation Functions" "$SHELL_RC" 2>/dev/null; then
        echo "Navigation functions already installed in $SHELL_RC"
        echo "To update, run: ./navigation_scripts.sh --uninstall && ./navigation_scripts.sh --install"
        exit 0
    fi

    echo "Installing navigation functions to $SHELL_RC..."
    echo "$NAVIGATION_FUNCTIONS" >> "$SHELL_RC"
    echo "✅ Navigation functions installed!"
    echo ""
    echo "To activate in current shell, run:"
    echo "  source $SHELL_RC"
    echo ""
    echo "Usage examples:"
    echo "  cd analysis/paper"
    echo "  s              # Jump to scripts/ from paper/"
    echo "  d              # Jump to data/ from scripts/"
    echo "  p              # Jump back to paper/"
    echo "  r              # Jump to project root"
    echo "  nav            # List all shortcuts"
}

# Function to uninstall navigation functions
uninstall_functions() {
    if ! grep -q "ZZCOLLAB Navigation Functions" "$SHELL_RC" 2>/dev/null; then
        echo "Navigation functions not found in $SHELL_RC"
        exit 0
    fi

    echo "Removing navigation functions from $SHELL_RC..."
    # Remove lines between markers (including markers)
    sed -i.bak '/# ZZCOLLAB Navigation Functions/,/# End ZZCOLLAB Navigation Functions/d' "$SHELL_RC"
    echo "✅ Navigation functions removed!"
    echo "Backup saved to: ${SHELL_RC}.bak"
    echo ""
    echo "To deactivate in current shell, run:"
    echo "  source $SHELL_RC"
}

# Main logic
case "$1" in
    --install)
        install_functions
        ;;
    --uninstall)
        uninstall_functions
        ;;
    *)
        echo "ZZCOLLAB Navigation Functions Setup"
        echo ""
        echo "This script installs shell functions for one-letter navigation"
        echo "that work from ANY subdirectory in your ZZCOLLAB project."
        echo ""
        echo "Usage:"
        echo "  ./navigation_scripts.sh --install    Install navigation functions"
        echo "  ./navigation_scripts.sh --uninstall  Remove navigation functions"
        echo ""
        echo "After installation, you can use:"
        echo "  r  → Jump to project root"
        echo "  d  → Jump to data/"
        echo "  s  → Jump to analysis/scripts/"
        echo "  p  → Jump to analysis/paper/"
        echo "  f  → Jump to analysis/figures/"
        echo "  nav → List all shortcuts"
        echo ""
        echo "Example workflow:"
        echo "  cd analysis/paper    # Working on paper"
        echo "  s                    # Jump to scripts to edit analysis"
        echo "  d                    # Jump to data to check raw data"
        echo "  p                    # Jump back to paper"
        ;;
esac
