package hep.lcio.implementation.sio;

import hep.io.sio.SIOInputStream;
import hep.io.sio.SIOOutputStream;
import hep.io.sio.SIORef;

import hep.lcio.event.TrackerHit;

import hep.lcio.implementation.event.ITrackerHit;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

/**
 *
 * @author Tony Johnson
 * @version $Id: SIOTrackerHit.java,v 1.12 2010/02/03 09:53:14 gaede Exp $
 */
class SIOTrackerHit extends ITrackerHit
{
   private List tempHits;
   SIOTrackerHit(SIOInputStream in, SIOEvent owner, int major, int minor) throws IOException
   {
      setParent(owner);
      type = in.readInt();
      for (int i = 0; i < 3; i++)
         position[i] = in.readDouble();
      for (int i = 0; i < 6; i++)
         covMatrix[i] = in.readFloat();
      dEdx = in.readFloat();
      time = in.readFloat();
      
      quality = 0 ; 
      if( SIOVersion.encode(major,minor) > SIOVersion.encode(1,11))
        quality = in.readInt();
      
      int nRawHits = 1 ;
      if( SIOVersion.encode(major,minor) > SIOVersion.encode(1,2))
      {
         nRawHits = in.readInt() ;
      }
      
      tempHits = new ArrayList( nRawHits ) ;
      for (int i = 0; i < nRawHits ; i++)
      {
         tempHits.add(in.readPntr());
      }
      rawHits = null;
      
      in.readPTag(this);
   }
   public List getRawHits()
   {
      if (rawHits == null && tempHits != null)
      {
         for (ListIterator i = tempHits.listIterator(); i.hasNext(); )
         {
            SIORef ref = (SIORef) i.next();
            i.set(ref.getObject());
         }
         rawHits = tempHits;
         tempHits = null;
      }
      return super.getRawHits();
   }   
   
   static void write(TrackerHit hit, SIOOutputStream out) throws IOException
   {
      if (hit instanceof SIOTrackerHit)
      {
         ((SIOTrackerHit) hit).write(out);
      }
      else
      {
         out.writeInt(hit.getType());
         double[] pos = hit.getPosition();
         for (int i = 0; i < 3; i++)
            out.writeDouble(pos[i]);
         
         float[] matrix = hit.getCovMatrix();
         for (int i = 0; i < 6; i++)
            out.writeFloat(matrix[i]);
         out.writeFloat(hit.getdEdx());
         out.writeFloat(hit.getTime());
         out.writeInt(hit.getQuality());
         
         List rawHits = hit.getRawHits() ;
         out.writeInt( rawHits.size()) ;
         for (int i = 0; i < rawHits.size() ; i++)
         {
            out.writePntr( rawHits.get(i) );
         }
         
         out.writePTag(hit);
      }
   }
   
   private void write(SIOOutputStream out) throws IOException
   {
      out.writeInt(type);
      for (int i = 0; i < 3; i++)
         out.writeDouble(position[i]);
      for (int i = 0; i < 6; i++)
         out.writeFloat(covMatrix[i]);
      out.writeFloat(dEdx);
      out.writeFloat(time);
      out.writeInt(quality);
      
      out.writeInt( rawHits.size()) ;
      for (int i = 0; i < rawHits.size() ; i++)
      {
         out.writePntr( rawHits.get(i) );
      }
      out.writePTag(this);
   }   
}
