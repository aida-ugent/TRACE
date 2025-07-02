import { useState, useRef, useLayoutEffect } from "react";

export function DefaultButton(props) {

  const { children, className, ...restProps } = props;

  return (
    <button
      type="button"
      className={`min-w-[100px] h-fit select-none justify-center items-center 
      rounded-md bg-slate-500 text-gray-50 mr-2 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
       border-0 drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-slate-800
       hover:opacity-75 hover:text-gray-50 focus:opacity-100 focus:shadow-none 
      focus:outline-none disabled:opacity-75 disabled:bg-gray-300 disabled:text-gray-500 ${className || ''}`}
      {...restProps}
    >
      {children}
    </button>
  )
}

export function SubmitButton(props) {

  const { children, className, ...restProps } = props;

  return (
    <button
      type="submit"
      className={`min-w-[100px] h-fit select-none justify-center items-center 
      rounded-md bg-slate-500 text-gray-50 mr-2 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
       border-0 drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-slate-800
       hover:opacity-75 hover:text-gray-50 focus:opacity-100 focus:shadow-none 
      focus:outline-none disabled:opacity-75 disabled:bg-gray-300 disabled:text-gray-500 ${className || ''}`}
      {...restProps}
    >
      {children}
    </button>
  )
}


export function AsyncButton(props) {

  const { onClick, children, className, ...restProps } = props;
  const [isLoading, setLoading] = useState(false);
  const [buttonWidth, setWidth] = useState("0");
  const ref = useRef();

  const test_dimensions = () => {
    if (ref.current) {
      setWidth((ref.current.offsetWidth + 1).toString());
    }
  }

  useLayoutEffect(() => {
    test_dimensions();
  }, []);

  const handleClick = () => {
    setLoading(true);
    onClick()
      .then(res => setLoading(false));
  }

  return (
    <>
      <button
        type="button"
        className={`min-w-[100px] h-fit select-none justify-center items-center 
        rounded-md bg-slate-500 text-gray-50 mr-2 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
         border-0 drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-slate-800
         hover:opacity-75 hover:text-gray-50 focus:opacity-100 focus:shadow-none 
        focus:outline-none disabled:opacity-75 disabled:bg-gray-300 disabled:text-gray-500 ${className || ''}`}
        onClick={handleClick}
        disabled={isLoading}
        ref={ref}
        {...restProps}
      >
        {isLoading ? "loading..." : children}
      </button>
    </>
  )

}

export function ChevronButton(props) {
  const { innerText, onChange, ...other } = props;
  const [visibility, setVisibility] = useState("hidden")

  const toggleVisibility = () => {
    if (visibility == "visible") {
      setVisibility("hidden")
    }
    else {
      setVisibility("visible")
    }
  }

  return (
    <>
      <DefaultButton onClick={toggleVisibility}>
        {<>
          <p className="select-none">{innerText}</p>
          {visibility == "visible" ? (<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-4 h-4 ml-2">
            <path fillRule="evenodd" d="M11.47 7.72a.75.75 0 011.06 0l7.5 7.5a.75.75 0 11-1.06 1.06L12 9.31l-6.97 6.97a.75.75 0 01-1.06-1.06l7.5-7.5z" clipRule="evenodd" />
          </svg>) : (<svg xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 24 24" fill="currentColor" className="w-4 h-4 ml-2">
            <path fillRule="evenodd" d="M12.53 16.28a.75.75 0 01-1.06 0l-7.5-7.5a.75.75 0 011.06-1.06L12 14.69l6.97-6.97a.75.75 0 111.06 1.06l-7.5 7.5z" clipRule="evenodd" />
          </svg>)}
        </>
        }
      </DefaultButton>
    </>
  )
}

export function ResetButton(props) {
  const { children, className, ...restProps } = props;
  
  return (
    <button
      type="button"
      className={`select-none inline-block rounded-full bg-white/80 p-2 text-sm font-medium leading-normal 
      text-white border-none hover:no-underline focus:opacity-100 focus:shadow-none focus:outline-none ${className || ''}`}
      {...restProps}
    >
      {children || (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6 fill-gray-500 hover:fill-gray-800">
          <path d="M6 3a3 3 0 0 0-3 3v1.5a.75.75 0 0 0 1.5 0V6A1.5 1.5 0 0 1 6 4.5h1.5a.75.75 0 0 0 0-1.5H6ZM16.5 3a.75.75 0 0 0 0 1.5H18A1.5 1.5 0 0 1 19.5 6v1.5a.75.75 0 0 0 1.5 0V6a3 3 0 0 0-3-3h-1.5ZM12 8.25a3.75 3.75 0 1 0 0 7.5 3.75 3.75 0 0 0 0-7.5ZM4.5 16.5a.75.75 0 0 0-1.5 0V18a3 3 0 0 0 3 3h1.5a.75.75 0 0 0 0-1.5H6A1.5 1.5 0 0 1 4.5 18v-1.5ZM21 16.5a.75.75 0 0 0-1.5 0V18a1.5 1.5 0 0 1-1.5 1.5h-1.5a.75.75 0 0 0 0 1.5H18a3 3 0 0 0 3-3v-1.5Z" />
        </svg>
      )}
    </button>
  )
}


export function ScreenshotButton(props) {
  const { children, className, ...restProps } = props;
  
  return (
    <button
      type="button"
      className={`select-none inline-block rounded-full bg-white/80 p-2 text-sm font-medium leading-normal 
      text-white border-none hover:no-underline focus:opacity-100 focus:shadow-none focus:outline-none ${className || ''}`}
      {...restProps}
    >
      {children || (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6 fill-gray-500 hover:fill-gray-800">
          <path d="M12 9a3.75 3.75 0 1 0 0 7.5A3.75 3.75 0 0 0 12 9Z" />
          <path fillRule="evenodd" d="M9.344 3.071a49.52 49.52 0 0 1 5.312 0c.967.052 1.83.585 2.332 1.39l.821 1.317c.24.383.645.643 1.11.71.386.054.77.113 1.152.177 1.432.239 2.429 1.493 2.429 2.909V18a3 3 0 0 1-3 3h-15a3 3 0 0 1-3-3V9.574c0-1.416.997-2.67 2.429-2.909.382-.064.766-.123 1.151-.178a1.56 1.56 0 0 0 1.11-.71l.822-1.315a2.942 2.942 0 0 1 2.332-1.39ZM6.75 12.75a5.25 5.25 0 1 1 10.5 0 5.25 5.25 0 0 1-10.5 0Zm12-1.5a.75.75 0 1 0 0-1.5.75.75 0 0 0 0 1.5Z" clipRule="evenodd" />
        </svg>
      )}
    </button>
  )
}



export function ZoomButton(props) {
  const { children, className, ...restProps } = props;
  
  return (
    <button
      type="button"
      className={`select-none inline-block rounded-full bg-gray-600 mx-2 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none ${className || ''}`}
      {...restProps}
    >
      {children || (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
          <path fillRule="evenodd" d="M10.5 3.75a6.75 6.75 0 100 13.5 6.75 6.75 0 000-13.5zM2.25 10.5a8.25 8.25 0 1114.59 5.28l4.69 4.69a.75.75 0 11-1.06 1.06l-4.69-4.69A8.25 8.25 0 012.25 10.5z" clipRule="evenodd" />
        </svg>
      )}
    </button>
  )
}

export function SettingsButton(props) {
  const { children, className, ...restProps } = props;
  
  return (
    <button
      type="button"
      className={`select-none inline-block bg-transparent mx-2 p-2 text-sm font-medium 
      leading-normal text-slate-400 border-none 
      hover:no-underline hover:text-slate-900
      focus:opacity-100 focus:shadow-none focus:outline-none ${className || ''}`}
      {...restProps}
    >
      {children || (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
          <path d="M18.75 12.75h1.5a.75.75 0 000-1.5h-1.5a.75.75 0 000 1.5zM12 6a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 0112 6zM12 18a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 0112 18zM3.75 6.75h1.5a.75.75 0 100-1.5h-1.5a.75.75 0 000 1.5zM5.25 18.75h-1.5a.75.75 0 010-1.5h1.5a.75.75 0 010 1.5zM3 12a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 013 12zM9 3.75a2.25 2.25 0 100 4.5 2.25 2.25 0 000-4.5zM12.75 12a2.25 2.25 0 114.5 0 2.25 2.25 0 01-4.5 0zM9 15.75a2.25 2.25 0 100 4.5 2.25 2.25 0 000-4.5z" />
        </svg>
      )}
    </button>
  )
}


export function ChevronRightButton(props) {
  const { children, className, ...restProps } = props;
  
  return (
    <button
      type="button"
      className={`select-none inline-block bg-transparent pt-3 px-1 text-sm font-medium 
      leading-normal text-slate-400 border-none 
      hover:no-underline hover:text-slate-900
      focus:opacity-100 focus:shadow-none focus:outline-none ${className || ''}`}
      {...restProps}
    >
      {children || (
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor" className="w-6 h-6">
          <path fillRule="evenodd" d="M15.28 9.47a.75.75 0 0 1 0 1.06l-4.25 4.25a.75.75 0 1 1-1.06-1.06L13.69 10 9.97 6.28a.75.75 0 0 1 1.06-1.06l4.25 4.25ZM6.03 5.22l4.25 4.25a.75.75 0 0 1 0 1.06l-4.25 4.25a.75.75 0 0 1-1.06-1.06L8.69 10 4.97 6.28a.75.75 0 0 1 1.06-1.06Z" clipRule="evenodd" />
        </svg>
      )}
    </button>
  )
}