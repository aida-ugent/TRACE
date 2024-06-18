import React, { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "./styles.css";

import Scatterplot from "./scatterplot";

window.addEventListener('error', e => {
  const { message, filename, lineno, colno, error } = e;

  if (message === "ResizeObserver loop completed with undelivered notifications.") {
    console.log("Window error: ResizeObserver loop completed with undelivered notifications.")
    const resizeObserverErr = document.getElementById('webpack-dev-server-client-overlay');
    if (resizeObserverErr) {
      resizeObserverErr.style.display = 'none';
    }
  }
});

const root = createRoot(document.getElementById("root"));
root.render(
  <StrictMode>
    <main className="h-screen w-screen flex flex-row items-top">
      <Scatterplot />
    </main>
  </StrictMode>
);