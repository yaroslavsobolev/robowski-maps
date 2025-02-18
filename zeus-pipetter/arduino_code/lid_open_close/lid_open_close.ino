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
