This code is supposed to be executed on computer
connected to CRAIC microspectrometer and get spectra
automatically via pyautogui interface to CRAIC software LabmdaFire

# Instructions
1. Put python window into designared area of the screen (consult desktop background).
2. Launch the Python script.
3. If the LambdaFire or Stage Control applications are not running yet, Python will detect it and do these things automatically:
   1. Launch the CRAIC Lambda Fire application (this window must have the smallest possible size both vertically and horizontally).
   2. Launch Craic Stage Control application.
   3. Click "Goto Abs" button on the Craic Stage Control and move the popped-up window so that it does not block any other windows.