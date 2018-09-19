package processing.test.copflatmap;

import com.modestmaps.InteractiveMap;

import ketai.ui.KetaiAlertDialog;
import processing.core.*;
import processing.data.*;

public class CoPFlatMap extends PApplet{ // implements View.OnClickListener {

    float[] x = new float[3600];
    float[] y = new float[3600];

    double GHA = Math.toRadians(106.076);
   // double dec = Math.toRadians(61.655);
    double dec = Math.toRadians(32.655);
 //   float GHA = 106.075f; //dubhe
 //   float dec = 61.655f;

    Table earthquakes, delta;

    int count;
    PImage world;

    double[]vv = new double[3]; //, vy[3], vyz[3];
    double[]vy = new double[3];
    double[]vyz = new double[3];
    double[]wpt = new double[3600];
    float[]WPT = new float[3600];

    double[][]My = new double[3][3];

    double[][]Mz = new double[3][3];

    String src =
            "http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_hour.csv";  // 1
    double longitude, latitude, altitude;
    //KetaiLocation location;                                   // 1
    float starX, starY;
    float accuracy;
    float myX, myY;
    public float[] Data;
    //   float alti = 90f;
    float alt = 90 - 12.8983f;
    float Alt;
    //   float alt1;
    PFont font;
    boolean Dubhe = false;
    boolean Procyon = false;
    float Dist = 0;
    float Dist1 = 0;

    float c1x;
    float c1y;
    float c2x,c2y;
    float r1,r2;
    float Lat1 = 39.24f;
    float Lat2;
    float Lon1 = -77.24f;
    float Lon2;
    float Alt1,Alt2;
    float Azi1,Azi2;

    InteractiveMap map;

    public void setup()
    {

        // location = new KetaiLocation(this);
  /*      try {

            earthquakes = loadTable(src, "header, csv");          // 2
        }
        catch
                (Exception x) {
            println("Failed to open online stream reverting to local data");
            earthquakes = loadTable("all_hour_2015-02-24.txt", "header, csv");      // 3
            earthquakes = loadTable("all_hour.csv", "header, csv");
        }
        count = earthquakes.getRowCount(); */
   //     println(count + " earthquakes found");

        orientation(LANDSCAPE);
        // world = loadImage("Worldmap.png");
        //  world = loadImage("Topo.jpg");
        //  world = loadImage("Time.png");
        world = loadImage("EquiSW.jpg");

        //EquiSW.jpg By Strebe - Own work, CC BY-SA 3.0, https://commons.wikimedia.org/w/index.php?curid=16115228
        //  world = loadImage("World_Time_Zones.png");
        font = loadFont("CourierNew36.vlw");


    }

    public void draw ()
    {
        background(127);
        image(world, 0, 0, width, height);                      // 4

        MainActivity.sendAttitudeRequest();

        Data = MainActivity.MyHandler.getData();
        x = new float[3600];
        y = new float[3600];

     //   Mz = Rz(Math.toRadians(360.0) - GHA, Mz);
        Mz = Rz(Math.toRadians(360.0) - Math.toRadians(Data[2]), Mz);



        My =  Ry(Math.toRadians(90.0) - dec, My);



        int w = 0;

        for( double L0 = -180.0; L0 <= 180.0; L0 += .1 )
        {
//float alt = map(mouseX,width-width,width,0,90);
            //  vv =  VectorSpherical2Cartesian(Math.toRadians(62.5),Math.toRadians(L0) );
            vv =  VectorSpherical2Cartesian(Math.toRadians(-Data[1]),Math.toRadians(L0) );


            vy =  MatrixVecProd( My, vv, vy );

            vyz =  MatrixVecProd( Mz, vy, vyz );



            wpt[w] = C2ELat( vyz[0], vyz[1], vyz[2]);
            wpt[w+1] = C2ELon( vyz[0], vyz[1], vyz[2]);

            WPT = toFloatArray(wpt);

            x[w] = map(WPT[w+1],radians(-180) ,radians(180),width, width - width);
            y[w] = map(WPT[w],radians(-90),(radians(90)),height,height - height);
            point(x[w],y[w]);
            noFill();
            stroke(0,0,0);
            strokeWeight(3);
            beginShape();
            curveVertex(x[w],y[w]);
            curveVertex(x[w],y[w]);

            endShape();
        }
        w++;


   /*     for (int row = 0; row < count; row++)
        {


            noFill();

            noStroke();
            fill(255, 127, 127);


            noStroke();
            myX = map(-77.24f, -180, 180, 0, width);
//   myY = map(39.24f, 85, -60, 0, height);
            myY = map(39.24f, 90, -90, 0, height);

            accuracy = 1f;

            noFill();
            stroke(150);
            strokeWeight(3f);
            //  sin(Ho) = sin(Lat) * sin(Dec) + cos(Lat) * cos(dec) * cos(LHA)
            //  LHA = GHA + Lon = 28.835
//28.835
            //.87885367288944529617935247451008 = sin(Ho)

            float GHA = 106.075f; //dubhe
            float dec = 61.655f;
            alt = 90 - 61.505f;
            float  GHA1 = map(GHA,0,360,0, width);
            //  float  dec1 = map (dec,85,-60,0,height);
            float  dec1 = map (dec,90,-90,0,height);
            //   float  alt1 = map(alt, 0, 90, 0, width);
            //   float  alt2 = map(alt, 0, 90, 0, height);
            ellipse(GHA1, dec1, dist(myX, myY, GHA1, dec1) * 2, dist(myX, myY, GHA1, dec1) * 2); //dubhe equal altitude
            ellipse(GHA1, dec1, dist(myX, myY, GHA1, dec1) * 2, dist(myX, myY, GHA1, dec1) * 2); //dubhe equal altitude
            alt = -Data[1];
            float alt1 = map((alt),90f,0f,0,width);
            float alt2 = map((alt),90f,0f,0,height);
            //  if (Dubhe == true) {
            //    ellipse(GHA1, dec1, dist(alt1, alt2, GHA1, dec1), dist(alt1, alt2, GHA1, dec1)); //live altitude
            //  }
            if (Dubhe == true) {
         /*       ellipse(GHA1, dec1, dist(alt1, alt2, GHA1, dec1), dist(alt1, alt2, GHA1, dec1)); //live altitude
                Dist1 = dist(alt1,alt2,GHA1,dec1);
                //    c1x = (int)GHA1;
                //    c1y = (int)dec1;
                r1= Dist1 / 2;
                float diag = dist(0,0,width,height);
                // float m1 = 60 * (90 - (-Data[1]));
                float m1 = 60 * (90 - (61.505f));
                float m2 = map(Dist1,0,diag,0,m1 ); */
            //    println("Nautical Miles  " + m1);

/*

            }
            if (Dubhe == false){
                ellipse(GHA1, dec1,Dist1,Dist1);


            }
            noStroke();
            fill(255, 0, 0); //red
            ellipse(GHA1, dec1, 10, 10); //dubhe GP
            //  ellipse(map(157.2217f, 0, 360, 0, width), map(5.1767f, 85, -60, 0, height), 10, 10);  //procyon GP
            ellipse(map(157.2217f, 0, 360, 0, width), map(5.1767f, 90, -90, 0, height), 10, 10);  //procyon GP
            noFill();
            GHA = 157.2217f; //procyon
            dec = 5.1767f;





            GHA1 = map(GHA,0,360,0, width);
            // dec1 = map (dec,85,-60,0,height);
            dec1 = map (dec,90,-90,0,height);
            starX = GHA1;
            starY = dec1;
            noFill();
            stroke(150);
            strokeWeight(3f);

            ellipse(GHA1, dec1, dist(myX, myY, GHA1, dec1) * 2, dist(myX, myY, GHA1, dec1) * 2); //procyon equal altitude
            if (Procyon == true) {
                ellipse(GHA1, dec1, dist(alt1, alt2, GHA1, dec1), dist(alt1, alt2, GHA1, dec1)); //live altitude
                Dist = dist(alt1,alt2,GHA1,dec1);
                //    c2x = (int)GHA1;
                //    c2y = (int)dec1;
                r2= Dist / 2;


            }
            if (Procyon == false){
                ellipse(GHA1,dec1,Dist,Dist);


            }

            if (Dubhe == false && Procyon == false) {
                intersections();
                noStroke();
                // Lon1 = map()
                fill(255, 0, 100);
                // ellipse((Lon1),Lat1,10,10);
                // ellipse(Lon2,Lat2,10,10);
                //  map.draw();
            }


            fill(0);

        }  */

        // Current Device location
        noStroke();
        float s = map(millis() % (100 * accuracy * 3.28f), 0, 100*accuracy * 3.28f, 0, 127); // 16
        fill(127, 255, 127);
        ellipse(myX, myY, 10, 10);
        fill(127, 255, 127, 127 - s);
        ellipse(myX, myY, s, s);



        textFont(font, 20);
        textAlign(LEFT, TOP);


        fill(0);
        noStroke();

        float rw = textWidth("Azimuth: " + Data[2] + "  " + "Altitude: " + (-Data[1]) ); //("Geodetic Position: " + -77.24f + "W " + 39.24f + "N   " );

        rect(width - 205 - rw, height - 65 - g.textSize, rw, g.textSize + textDescent());
        fill(255, 255, 0);
        textAlign(RIGHT, BOTTOM);

        text("Azimuth: " + Data[2] + "  " + "Altitude: " + (-Data[1]), width - 205, height - 65);
        // println("rw: " + rw);

        fill(0);
        noStroke();
        rw = textWidth("Lat: " + Lat1 + "  " + "Lon: " + Lon1 );
        rect(width - 205 - rw, height - 85 - g.textSize, rw, g.textSize + textDescent());
        rect(width - 205 - rw, height - 105 - g.textSize, rw, g.textSize + textDescent());
        fill(255, 255, 0);
        textAlign(RIGHT, BOTTOM);

        text("Lat: " + Lat1 + "  " + "Lon: " + Lon1 , width-205, height-85);
        text("Lat: " + Lat2 + "  " + "Lon: " + Lon2 , width-205, height-105);


    }

    public void mousePressed()
    {
      /*  if (mouseY < (starY + 100))
        {
            if (mouseX < width/3)
                KetaiKeyboard.toggle(this);
          //  else if (mouseX > width/3 && mouseX < width-(width/3))
             //   KetaiAlertDialog.popup(this, "Pop Up!", "this is a popup message box");
            else if (mouseX > starX && mouseX < width-(starX) )
                KetaiAlertDialog.popup(this, "Altitude", "Point at star and click OK");

          //  else if (mouseX > (starX -100) && mouseX < width-(starX - 100) )
           //     KetaiAlertDialog.popup(this, "Azimuth", "Point at star and click OK");
             //   vibe.vibrate(1000);
        } */
        if(mouseX < (map(106.075f,0,360,0, width)+ 10) && mouseX > (map(106.075f, 0, 360, 0, width) - 10)) {

            KetaiAlertDialog.popup(this, "Altitude", "Point at star and click OK");
            Dubhe = true;
            Procyon = false;


        }
        if(mouseX < (map(157.2217f,0,360,0, width)+ 10) && mouseX > (map(157.2217f, 0, 360, 0, width) - 10)) {

            KetaiAlertDialog.popup(this, "Altitude", "Point at star and click OK");
            Dubhe = false;
            Procyon = true;


        }

        if (mouseX > width/2){
            Dubhe = false;
            Procyon = false;
        }

    }

    void intersections() {
        //    float d = dist(c1x, c1y, c2x, c2y); // distance between centers

        c1x =   map( 106.075f,0,360,0, width);
        // c1y =  map (61.655f,85,-60,0,height);
        c1y =  map (61.655f,90,-90,0,height);
        c2x =   map(157.2217f,0,360,0, width);

        // c2y =  map (5.1767f,85,-60,0,height);
        c2y =  map (5.1767f,90,-90,0,height);
        float d = dist(c1x, c1y, c2x, c2y); // distance between centers
        //  float d = dist(map( 106.075f,0,360,0, width),map (61.655f,85,-60,0,height),map(157.2217f,0,360,0, width),map (5.1767f,85,-60,0,height));
        float base, h; // auxiliary distances

        //  p, middle point between q1 and q2
        // q1 dn q2 intersection points
        float px, py, q1x, q1y, q2x, q2y;


        //   else{ // intersect in two points

        base=(r1*r1-r2*r2+d*d)/(2*d);
        h=sqrt(r1*r1-base*base);

        px=c1x+base*(c2x-c1x)/d;
        py=c1y+base*(c2y-c1y)/d;


        q1x=px+h*(c2y-c1y)/d;
        q1y=py-h*(c2x-c1x)/d;
        q2x=px-h*(c2y-c1y)/d;
        q2y = py +h*(c2x-c1x) / d;

        noStroke();
        fill(0, 0, 255);
        //     ellipse(px, py, 10, 10);
        ellipse(q1x,q1y,10,10);
        ellipse(q2x,q2y,10,10);
        // Lat1 = map(q2y,0,height,85,-60);
        Lat1 = map(q2y,0,height,90,-90);
        Lon1 = map(q2x,0,width,-180,180);
        // Lat1 = map(q2y,0,height,85,-60);
        Lat1 = map(q2y,0,height,90,-90);
        Lon1 = map(q2x,0,width,-180,180);
        Lat2 = map(q1y,0,height,90,-90);
        Lon2 = map(q1x,0,width,-180,180);
        //   myX = map(-77.24f, -180, 180, 0, width);
        //   myY = map(39.24f, 85, -60, 0, height);

        //       println("intersections: Q1=("+ q1x+ ", "+ q1y+ ") and  Q2=("+q2x+ ", "+ q2y+")");
        //    println("Position1: Lat = " + ( (Lat1)) + " " + "Lon = " + Lon1);
        //  }


        
    }

    //Funcition to convert double[] to float[]
    float[] toFloatArray(double[] arr) {
        if (arr == null) return null;
        int n = arr.length;
        float[] ret = new float[n];
        for (int i = 0; i < n; i++) {
            ret[i] = (float)arr[i];
        }
        return ret;
    }
// end of function to convert double[] to float[]

  
    
    
/*
File: vector.cpp
Cálculo vectorial
Resultado test: OK
This file contains proprietary information of Andrés Ruiz Gonzalez
Andrés Ruiz. San Sebastian - Donostia. Gipuzkoa
Copyright (c) 1999 - 2007
*/    
    
    
    
    
    
    
    
    
    
    double[] VectorSpherical2Cartesian(double B, double L){

        double v[] = new double[3];
        v[0] = Math.cos(B) * Math.cos(L);
        v[1] = Math.cos(B) * Math.sin(L);
        v[2] = Math.sin(B);
//   println(B);
//   println(L);
        return(v);

    }

    public double C2ELat( double x, double y, double z )
    {
        double[]res = new double[3];
        res[0] = Math.sqrt( x*x+y*y+z*z);  //R
//*B = ASIN(z/(*R));
        res[1] = Math.atan2( z, Math.sqrt(x*x+y*y) ); //B
        res[2] = Math.atan2( y, x ); //L

//println("R:  " + (res[0]) + "  B: " + Math.toDegrees(res[1]) + "  L: " + Math.toDegrees(res[2]));

        return (res[1]);
//println(R);
    }

    public double C2ELon( double x, double y, double z )
    {
        double[]res = new double[3];
        res[0] = Math.sqrt( x*x+y*y+z*z);  //R
//*B = ASIN(z/(*R));
        res[1] = Math.atan2( z, Math.sqrt(x*x+y*y) ); //B
        res[2] = Math.atan2( y, x ); //L

//println("R:  " + (res[0]) + "  B: " + Math.toDegrees(res[1]) + "  L: " + Math.toDegrees(res[2]));

        return (res[2]);
//println(R);
    }

//public double[] E2C( double B, double L, double R, double x, double y, double z )

    public double[] E2C( double B, double L, double R )
    {
        double[]res = new double[3];

        res[0] = R*Math.cos((B))*Math.cos((L));
        res[1] = R*Math.cos((B))*Math.sin((L));
        res[2] = R*Math.sin((B));

// println(res[0] + " " + res[1] + " " + res[0]);

        return(res);
    }

    public double[][] Rx( double a, double[][] M ){

        M[0][0] = 1.0;
        M[1][0] = 0.0;
        M[2][0] = 0.0;
        M[0][1] = 0.0;
        M[1][1] = Math.cos(a); //Math.cos(Math.toRadians(a));
        M[2][1] = Math.sin(a); //Math.sin(Math.toRadians(a));
        M[0][2] = 0.0;
        M[1][2] = -Math.sin(a); //-Math.sin(Math.toRadians(a));
        M[2][2] = Math.cos(a); //Math.cos(Math.toRadians(a));

        return(M);
    }

    public double[][] Ry( double a, double[][] M ){

        M[0][0] = Math.cos(a);
        M[1][0] = 0.0;
        M[2][0] = -Math.sin(a);
        M[0][1] = 0.0;
        M[1][1] = 1.0;
        M[2][1] = 0.0;
        M[0][2] = Math.sin(a);
        M[1][2] = 0.0;
        M[2][2] = Math.cos(a);

        return(M);
    }

    public double[][] Rz( double a, double[][] M ){

        M[0][0] = Math.cos(a); //Math.cos(a);
        M[1][0] = Math.sin(a);
        M[2][0] = 0.0;
        M[0][1] = -Math.sin(a);
        M[1][1] = Math.cos(a);
        M[2][1] = 0.0;
        M[0][2] = 0.0;
        M[1][2] = 0.0;
        M[2][2] = 1.0;

        return(M);
    }

    public double[] MatrixVecProd( double[][] A, double[] v, double[] res ) {

        int i,j;
        int n = 3;

        for( i=0; i<n; i++ ) {
            res[i] = 0.0;
            for( j=0; j<n; j++ ) {
                res[i] += A[i][j]*v[j];

            }
        }

        return (res);
    }


    void onLocationEvent(double _latitude, double _longitude, double _altitude)
    {
        longitude = _longitude;
        latitude = _latitude;
        altitude = _altitude;
        println("lat/lon/alt: " + latitude + "/" + longitude + "/" + altitude);
    }
    static public void main(String[] passedArgs) {
        String[] appletArgs = new String[] { "CoPFlatMap" };
        if (passedArgs != null) {
            PApplet.main(concat(appletArgs, passedArgs));
        } else {
            PApplet.main(appletArgs);
        }
    }


}
