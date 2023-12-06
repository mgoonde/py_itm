This is the source of itm library `libitm.so`.

Each template structure is an istance of `t_template` type, containing details of the structure.

All templates which are added form a linked list, which is a member of the `t_itm_ptr` type.

Requirements:
 - compiled `IRA` library, file: `libira.so` is needed.

To compile this directory:
```
make
```
This should produce `../lib/libitm.so`.
