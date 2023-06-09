import java.awt.Color;
import java.awt.Graphics;
 
import javax.swing.JFrame;
import javax.swing.JPanel;
 
public class leniamain{
    public static void main(String[] args) {
        GameWindow gw = new GameWindow("Lenia",600,600);
        
        lenia life=new lenia(new double[]{1.0},0.15,0.017,26,5,600,600,1);
        
        /*
        smoothlife life=new smoothlife(6,18,0.257,0.365,
        0.336,0.549,0.147,0.028,1.0,600,600,1);
        */
        //life.rand_cir_arr(500);
        //life.rand_state();
        for(int i=0;i<55;i++){
            life.make_orbium_at((i/10)*100, (i%10)*100);
        }
        DrawCanvas cd=new DrawCanvas();
        gw.add(cd);
        cd.set_life(life);
        gw.setVisible(true);
        gw.startGameLoop();
    }
}

//スレッドとして利用したい->Runnableインタフェースが必要
class GameWindow extends JFrame implements Runnable{
    private Thread th = null;
    public GameWindow(String title, int width, int height) {
        super(title);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setSize(width,height);
        setLocationRelativeTo(null);
        setResizable(false);
    }
 
    //ゲームループの開始メソッド
    public synchronized void startGameLoop(){
        if ( th == null ) {
            th = new Thread(this);
            th.start();
        }
    }
    //ゲームループの終了メソッド
    public synchronized void stopGameLoop(){
        if ( th != null ) {
            th = null;
        }
    }
    public void run(){
        //ゲームループ（定期的に再描画を実行）
        while(th != null){
            try{
                Thread.sleep(25);
                repaint();
            }catch(InterruptedException e){
                e.printStackTrace();
            }
        }
    }
}
 
class DrawCanvas extends JPanel{
    lenia life;
    void set_life(lenia L){
        life=L;
    }
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        g.setColor(Color.BLACK);
        g.fillRect(0,0,600,600);
        life.disp_life(g);
        //life.plot_kernel(g, life.ker,life.R);
        life.calc_update();
        life.state_update();
    }
}