import '@/styles/globals.css'
import { Inter } from 'next/font/google'
import { Roboto } from "next/font/google";

const roboto = Roboto({ weight: "400", subsets: ["latin"] });
const inter = Inter({ subsets: ['latin'] })

export const metadata = {
  title: 'Trace - Exploring 2D Embeddings',
  description: 'Interactive visualization and quality analysis',
}

export default function RootLayout({ children }) {
  return (
    <html lang="en">
      <body className={roboto.className}>{children}</body>
    </html>
  )
}