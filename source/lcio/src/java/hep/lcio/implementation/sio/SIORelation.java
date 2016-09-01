package hep.lcio.implementation.sio;

import hep.io.sio.SIOInputStream;
import hep.io.sio.SIOOutputStream;
import hep.io.sio.SIORef;
import hep.lcio.event.LCObject;
import hep.lcio.event.LCRelation;
import hep.lcio.implementation.event.ILCRelation;
import java.io.IOException;

/**
 *
 * @author Tony Johnson
 * @version $Id: SIORelation.java,v 1.4 2007/11/07 20:46:23 jeremy Exp $
 */
class SIORelation extends ILCRelation
{
   private SIORef f;
   private SIORef t;
   
   SIORelation(SIOInputStream in, SIOEvent owner, int flag, int major, int minor) throws IOException
   {
      super();
      f = in.readPntr();
      t = in .readPntr();
      float weight = 1;
      if ((flag & 1<<31) != 0) weight = in.readFloat();
   }
   static void write(LCRelation relation, SIOOutputStream out, int flag) throws IOException
   {
      out.writePntr(relation.getFrom());
      out.writePntr(relation.getTo());
      if ((flag & 1<<31) != 0) out.writeFloat(relation.getWeight());
   }
   
   public LCObject getTo()
   {
      if (from == null) from = (LCObject) f.getObject();
      return from;
   }
   
   public LCObject getFrom()
   {
      if (to == null) to = (LCObject) t.getObject();
      return to;
   }
   
}
