
import smtplib, ssl, threading, os, sys, getopt

stamps = {"file1": 0}



def check_file(path,pwd):
    threading.Timer(1, check_file,args=[path,pwd]).start()    
    date=os.path.getmtime(path)
    old_date=stamps["file1"]
    if (date!=old_date):
        stamps["file1"]=date
        alert=path + " changed " + str(date-old_date) + " sec ago."
        send_mail(alert,pwd)
        print(alert)


        
        
    


def send_mail(alert,pwd):
    sender_email = "flori_n@yahoo.de"
    receiver_email = "florian.niedermann@gmail.com"
    message = "Subject: Code alert: "+alert+"\n\n"
    

    port = 465  # For SSL
    

    # Create a secure SSL context
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL("smtp.mail.yahoo.com", port, context=context) as server:
        server.login("flori_n@yahoo.de", pwd)
        server.sendmail(sender_email, receiver_email, message)
    


def main(argv):
   inputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile="])
   except getopt.GetoptError:
      print('alert.py -i <inputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('alert.py -i <inputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
   print ('Input file is: ', inputfile)
   pwd=input("Password:")
   stamps["file1"]=os.path.getmtime(inputfile)
   check_file(inputfile,pwd)
   



if __name__ == "__main__":
   main(sys.argv[1:])



