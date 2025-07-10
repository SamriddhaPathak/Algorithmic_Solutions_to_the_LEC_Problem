# Error
# print("Hello World ) -> Syntax error

# num1 = int(input("Enter a number:"))
# num2 = int(input("Enter another number:"))

# sum = num1 - num2 # This represents logical error
# print(f"sum = {sum}")

#Exceptions 
# x = int(input("Enter a number: "))
# print(x) - > to demonstrate ValueError

# try:
#     x = int(input("Enter a number: "))

# except ValueError:
#     print("Please Enter an Integer!!")
    
# print(f"number = {x}")

# try:
#     x = int(input("Enter a number: "))

# except ValueError:
#     print("Please Enter an Integer!!")

# else:
#     print(f"number = {x}")

# while True:
#     try:
#         x = int(input("Enter a number: "))
#     except ValueError:
#         pass
#     else:
#         break
# print(f"The entered number is: {x}")

'''
try:



except:
.......
.......
.......

else: 

finally: 


'''
# Division By Zero Error
# Task -> Program -> takes a number -> divides 10 by the entered number.
# while True:
#     try:
#         x = int(input("Enter a number: "))
#         result = 10/x
        
#     except ValueError:
#         print(ValueError)  
        
#     except ZeroDivisionError:
#         print(ZeroDivisionError)   
    
#     else: 
#         print(f"Result = {result}")
#         break

# def main():
#     x = get_int()
#     print(f"number = {x}")
    
    
# def get_int():
#     while True:
#         try:
#             x = int(input("Enter a number: "))
            
#         except ValueError:
#             pass
        
#         else: 
#             return x
        
# main()

#program = take 3 numbers as input -> Sum, product and operation = (num1+num2)/num3 by using functional programming
# while True:
#     try:
#         num1 = int(input("Enter a number: "))
#         num2 = int(input("Enter another number: "))
#         num3 = int(input("Enter another number: "))
        
#         operation = (num1+num2)/num3
        
#     except (ValueError, ZeroDivisionError):
#         pass

#     else:
#         sum = num1 + num2 + num3
#         product = num1 * num2 * num3
#         print(f"sum = {sum}")
#         print(f"Product = {product}")
#         print(f"Operation = {operation}")
#         break
          
def main(num1,num2,num3):
    print(f"The sum is: {num1 + num2 + num3}")
    print(f"The product is: {num1 * num2 * num3}")
    try:
        print(f"The operation (num1 + num2)/num3 is {(num1 + num2)/num3}")
    except ZeroDivisionError:
        print(f"num3 = {num3} thus Division Invalid!!")
    
def cal():
    while True:
        try:
            num1 = int(input("Enter first number: "))
            num2 = int(input("Enter second number: "))
            num3 = int(input("Enter third number: "))
            
        except ValueError:
            print("Please enter an integer!!")
        else:
            return main(num1,num2,num3)
cal()