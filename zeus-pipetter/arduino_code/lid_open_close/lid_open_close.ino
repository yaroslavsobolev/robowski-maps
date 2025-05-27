/**
 * Arduino Servo & Digital IO Controller
 * -------------------------------------
 * This sketch interfaces with a computer over a serial connection (9600 baud)
 * and listens for string-based commands to control:
 *  - A servo motor (attached to pin 9) mounted on the lid of a Nanodrop spectrophotometer
 *  - Three digital output pins (2, 3, 4)
 *
 * Functionality:
 * --------------
 * - "1"  → Opens the servo (sweeps from 30° to 140°)
 * - "0"  → Closes the servo (sweeps from 140° to 30°)
 *
 * - "21" → Sets digital pin 2 HIGH
 * - "20" → Sets digital pin 2 LOW
 * - "31" → Sets digital pin 3 HIGH
 * - "30" → Sets digital pin 3 LOW
 * - "41" → Sets digital pin 4 HIGH
 * - "40" → Sets digital pin 4 LOW
 *
 * Servo Behavior:
 * ---------------
 * - Servo is initialized to 30° on setup.
 * - Servo movement uses linear sweep with delays:
 *    - Opening: fast sweep (10ms step delay)
 *    - Closing: slower sweep (25ms step delay), for avoiding smashing the lid
 *
 * Setup:
 * ------
 * - Initializes Serial at 9600 baud
 * - Configures pins 2, 3, 4 as OUTPUT and sets them HIGH by default
 * - Attaches the servo to pin 9
 *
 * Notes:
 * ------
 * - Works well with Python scripts sending string commands
 * - Useful for controlling mechanical setups (valves, lights, etc.)
 *   where servo and digital output are triggered via software commands.
 */


#include <Servo.h>
Servo servo;
int angle = 0;

void setup() {
  Serial.begin(9600);
  Serial.setTimeout(25);
  
  pinMode(2, OUTPUT);
  digitalWrite(2, HIGH);
  pinMode(3, OUTPUT); 
  digitalWrite(3, HIGH);
  pinMode(4, OUTPUT); 
  digitalWrite(4, HIGH);

  servo.attach(9);
  servo.write(30); 
   
}

void close();
void open();

void loop()
{
while (Serial.available() == 0) {}
String user_input = Serial.readString();

user_input.trim();
//Serial.println(user_input);
if(user_input == "1"){
  //Serial.println("1");
  open();
}
else if(user_input == "0"){
  //Serial.println("0");
  close();
}
else if(user_input == "21"){
  digitalWrite(2, HIGH);
  }
else if(user_input == "20"){
  digitalWrite(2, LOW);
  }
 else if(user_input == "31"){
  digitalWrite(3, HIGH);
  }
else if(user_input == "30"){
  digitalWrite(3, LOW);
  }
else if(user_input == "41"){
  //Serial.println("41");
  digitalWrite(4, HIGH);
  }
else if(user_input == "40"){
  //Serial.println("40");
  digitalWrite(4, LOW);
  }
}

void close()
{
for(angle = 140; angle >30; angle-=1)    
{                                
  servo.write(angle);           
  delay(25); 
}
}

void open()
{
for(angle = 30; angle < 140; angle+=1)    
{                                
  servo.write(angle);           
  delay(10); 
}
}
