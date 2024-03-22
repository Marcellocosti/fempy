import os
# If you want to import alice3 sw add `export FEMPY_ALICE3=true` to your .bashrc
if os.getenv("FEMPY_ALICE3", 'false').lower() in ('true', '1'):
    from . import alice3
